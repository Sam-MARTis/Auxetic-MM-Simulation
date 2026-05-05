from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import ezdxf
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull, Delaunay


DXF_FILE_PATH = Path("dxf_files/pattern_1x1.dxf")
OUTPUT_IMAGE_PATH = Path("testing_meshing_preview.png")
POINT_SPACING = 0.08
FIGURE_WIDTH = 12.0
FIGURE_HEIGHT = 12.0
BOUNDARY_COLOR = "#111827"
BOUNDARY_WIDTH = 2.2
INTERIOR_LINE_COLOR = "#9ca3af"
INTERIOR_LINE_WIDTH = 0.8
TRIANGLE_EDGE_COLOR = "#60a5fa"
TRIANGLE_EDGE_WIDTH = 0.45
NODE_COLOR = "#ef4444"
NODE_SIZE = 16.0
BACKGROUND_COLOR = "white"


@dataclass(frozen=True, slots=True)
class MeshNode:
    index: int
    x: float
    y: float


def load_dxf_modelspace(*, file_path: Path) -> tuple[Any, Any]:
    document = ezdxf.readfile(str(file_path))
    return document.modelspace(), document


def entity_to_segment(*, entity: Any) -> tuple[np.ndarray, np.ndarray]:
    start = entity.dxf.start
    end = entity.dxf.end
    return (
        np.array([float(start.x), float(start.y)], dtype=float),
        np.array([float(end.x), float(end.y)], dtype=float),
    )


def extract_boundary_segments(*, modelspace: Any) -> list[tuple[np.ndarray, np.ndarray]]:
    boundary_segments: list[tuple[np.ndarray, np.ndarray]] = []
    for line in modelspace.query("LINE"):
        boundary_segments.append(entity_to_segment(entity=line))
    return boundary_segments


def extract_all_segments(*, modelspace: Any) -> list[tuple[np.ndarray, np.ndarray]]:
    segments: list[tuple[np.ndarray, np.ndarray]] = []

    for insert in modelspace.query("INSERT"):
        block = insert.block()
        matrix = insert.matrix44()
        for line in block.query("LINE"):
            start = matrix.transform(line.dxf.start)
            end = matrix.transform(line.dxf.end)
            segments.append(
                (
                    np.array([float(start.x), float(start.y)], dtype=float),
                    np.array([float(end.x), float(end.y)], dtype=float),
                )
            )

    segments.extend(extract_boundary_segments(modelspace=modelspace))
    return segments


def ordered_boundary_vertices(*, boundary_segments: Sequence[tuple[np.ndarray, np.ndarray]]) -> list[np.ndarray]:
    boundary_points: list[np.ndarray] = []
    for start, end in boundary_segments:
        boundary_points.append(start)
        boundary_points.append(end)

    if not boundary_points:
        raise ValueError("No boundary segments were found in the DXF file.")

    unique_points = np.unique(np.round(np.asarray(boundary_points, dtype=float), decimals=12), axis=0)
    if unique_points.shape[0] < 3:
        raise ValueError("Could not recover a closed boundary loop from the DXF linework.")

    hull = ConvexHull(unique_points)
    return [unique_points[index] for index in hull.vertices]


def point_in_polygon(*, point: np.ndarray, polygon: Sequence[np.ndarray]) -> bool:
    x_coord = float(point[0])
    y_coord = float(point[1])
    polygon_points = np.asarray(polygon, dtype=float)
    inside = False
    previous_index = len(polygon_points) - 1

    for current_index in range(len(polygon_points)):
        current_x = float(polygon_points[current_index][0])
        current_y = float(polygon_points[current_index][1])
        previous_x = float(polygon_points[previous_index][0])
        previous_y = float(polygon_points[previous_index][1])

        intersects = (current_y > y_coord) != (previous_y > y_coord)
        if intersects:
            intersection_x = (previous_x - current_x) * (y_coord - current_y) / ((previous_y - current_y) + 1.0e-30) + current_x
            if x_coord < intersection_x:
                inside = not inside

        previous_index = current_index

    return inside


def build_interior_points(*, boundary_vertices: Sequence[np.ndarray], spacing: float) -> np.ndarray:
    boundary_array = np.asarray(boundary_vertices, dtype=float)
    min_x = float(np.min(boundary_array[:, 0]))
    max_x = float(np.max(boundary_array[:, 0]))
    min_y = float(np.min(boundary_array[:, 1]))
    max_y = float(np.max(boundary_array[:, 1]))

    x_values = np.arange(min_x, max_x + spacing, spacing, dtype=float)
    y_values = np.arange(min_y, max_y + spacing, spacing, dtype=float)
    grid_x, grid_y = np.meshgrid(x_values, y_values)
    grid_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))

    inside_mask = np.array(
        [point_in_polygon(point=point, polygon=boundary_array) for point in grid_points],
        dtype=bool,
    )
    return grid_points[inside_mask]


def combine_unique_points(*, boundary_vertices: Sequence[np.ndarray], interior_points: np.ndarray) -> np.ndarray:
    combined_points = np.vstack((np.asarray(boundary_vertices, dtype=float), np.asarray(interior_points, dtype=float)))
    rounded_points = np.round(combined_points, decimals=12)
    return np.unique(rounded_points, axis=0)


def build_mesh_nodes(*, points: np.ndarray) -> list[MeshNode]:
    nodes: list[MeshNode] = []
    for index, point in enumerate(points.tolist()):
        nodes.append(MeshNode(index=index, x=float(point[0]), y=float(point[1])))
    return nodes


def triangulate_nodes(*, nodes: Sequence[MeshNode]) -> Delaunay:
    point_array = np.array([[node.x, node.y] for node in nodes], dtype=float)
    return Delaunay(point_array)


def triangle_is_inside_boundary(*, triangle_points: np.ndarray, boundary_vertices: Sequence[np.ndarray]) -> bool:
    centroid = np.mean(triangle_points, axis=0)
    return point_in_polygon(point=centroid, polygon=boundary_vertices)


def filter_triangles(*, triangulation: Delaunay, nodes: Sequence[MeshNode], boundary_vertices: Sequence[np.ndarray]) -> np.ndarray:
    point_array = np.array([[node.x, node.y] for node in nodes], dtype=float)
    kept_triangles: list[np.ndarray] = []
    for triangle in triangulation.simplices:
        triangle_points = point_array[triangle]
        if triangle_is_inside_boundary(triangle_points=triangle_points, boundary_vertices=boundary_vertices):
            kept_triangles.append(triangle)
    return np.asarray(kept_triangles, dtype=int)


def plot_mesh(*, boundary_vertices: Sequence[np.ndarray], nodes: Sequence[MeshNode], triangles: np.ndarray, all_segments: Sequence[tuple[np.ndarray, np.ndarray]]) -> None:
    point_array = np.array([[node.x, node.y] for node in nodes], dtype=float)
    figure, axes = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    figure.patch.set_facecolor(BACKGROUND_COLOR)
    axes.set_facecolor(BACKGROUND_COLOR)

    for start, end in all_segments:
        axes.plot(
            [float(start[0]), float(end[0])],
            [float(start[1]), float(end[1])],
            color=INTERIOR_LINE_COLOR,
            linewidth=INTERIOR_LINE_WIDTH,
            alpha=0.85,
            zorder=1,
        )

    for triangle in triangles:
        triangle_points = point_array[triangle]
        closed_triangle = np.vstack((triangle_points, triangle_points[0]))
        axes.plot(
            closed_triangle[:, 0],
            closed_triangle[:, 1],
            color=TRIANGLE_EDGE_COLOR,
            linewidth=TRIANGLE_EDGE_WIDTH,
            alpha=0.95,
            zorder=2,
        )

    boundary_array = np.asarray(boundary_vertices, dtype=float)
    closed_boundary = np.vstack((boundary_array, boundary_array[0]))
    axes.plot(
        closed_boundary[:, 0],
        closed_boundary[:, 1],
        color=BOUNDARY_COLOR,
        linewidth=BOUNDARY_WIDTH,
        zorder=3,
    )

    axes.scatter(
        point_array[:, 0],
        point_array[:, 1],
        s=NODE_SIZE,
        color=NODE_COLOR,
        zorder=4,
    )
    axes.set_aspect("equal", adjustable="box")
    axes.set_title("DXF mesh preview with nodes", fontsize=14)
    axes.set_xlabel("X")
    axes.set_ylabel("Y")
    axes.grid(False)
    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE_PATH, dpi=220, bbox_inches="tight")
    plt.show()


def main() -> None:
    modelspace, _document = load_dxf_modelspace(file_path=DXF_FILE_PATH)
    boundary_segments = extract_boundary_segments(modelspace=modelspace)
    all_segments = extract_all_segments(modelspace=modelspace)
    boundary_vertices = ordered_boundary_vertices(boundary_segments=boundary_segments)
    interior_points = build_interior_points(boundary_vertices=boundary_vertices, spacing=POINT_SPACING)
    mesh_points = combine_unique_points(boundary_vertices=boundary_vertices, interior_points=interior_points)
    nodes = build_mesh_nodes(points=mesh_points)
    triangulation = triangulate_nodes(nodes=nodes)
    triangles = filter_triangles(triangulation=triangulation, nodes=nodes, boundary_vertices=boundary_vertices)
    print(f"Boundary vertices: {len(boundary_vertices)}")
    print(f"Interior points: {len(interior_points)}")
    print(f"Mesh nodes: {len(nodes)}")
    print(f"Kept triangles: {len(triangles)}")
    plot_mesh(boundary_vertices=boundary_vertices, nodes=nodes, triangles=triangles, all_segments=all_segments)


if __name__ == "__main__":
    main()
