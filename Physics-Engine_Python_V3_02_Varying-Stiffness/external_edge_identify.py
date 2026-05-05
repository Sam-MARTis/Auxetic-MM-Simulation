from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import ezdxf
import matplotlib.pyplot as plt
import numpy as np


DXF_FILE_PATH = Path("dxf_files/pattern_1x1.dxf")
OUTPUT_IMAGE_PATH = Path("external_edge_identify_preview.png")
FIGURE_SIZE = (10.0, 10.0)
EXTERNAL_EDGE_COLOUR = "#dc2626"
EXTERNAL_EDGE_COLOURS = ["#dc2626", "#ef4444", "#f87171", "#fca5a5", "#fee2e2", "#b91c1c", "#991b1b", "#7f1d1d", "#6b1212", "#4c0e0e"]
INTERNAL_EDGE_COLOUR = "#16a34a"
EXTERNAL_EDGE_WIDTH = 2.6
INTERNAL_EDGE_WIDTH = 1.4


def parse_arguments() -> Path:
	parser = argparse.ArgumentParser(description="Display external edges in red and internal edges in green for a DXF file.")
	parser.add_argument("dxf_file_path", nargs="?", type=Path, default=DXF_FILE_PATH, help="Path to a DXF file from the provided set.")
	arguments = parser.parse_args()
	return arguments.dxf_file_path


def load_modelspace(file_path: Path) -> Any:
	document = ezdxf.readfile(str(file_path))
	return document.modelspace()


def point_xy(point: Any) -> np.ndarray:
	return np.array([float(point.x), float(point.y)], dtype=float)


def collect_edges(modelspace: Any) -> list[tuple[np.ndarray, np.ndarray, bool]]:
	edges: list[tuple[np.ndarray, np.ndarray, bool]] = []

	for line in modelspace.query("LINE"):
		edges.append((point_xy(point=line.dxf.start), point_xy(point=line.dxf.end), True))

	for insert in modelspace.query("INSERT"):
		block = insert.block()
		matrix = insert.matrix44()
		for line in block.query("LINE"):
			start = matrix.transform(line.dxf.start)
			end = matrix.transform(line.dxf.end)
			edges.append((point_xy(point=start), point_xy(point=end), False))

	return edges


def plot_edges(edges: list[tuple[np.ndarray, np.ndarray, bool]]) -> None:
	figure, axes = plt.subplots(figsize=FIGURE_SIZE)
	figure.patch.set_facecolor("white")
	axes.set_facecolor("white")

	external_label_used = False
	internal_label_used = False

	for start, end, is_external in edges:
		if is_external:
			color = EXTERNAL_EDGE_COLOURS[np.random.randint(0, len(EXTERNAL_EDGE_COLOURS))]
			width = EXTERNAL_EDGE_WIDTH
			label = "External edges" if not external_label_used else None
			external_label_used = True
		else:
			color = INTERNAL_EDGE_COLOUR
			width = INTERNAL_EDGE_WIDTH
			label = "Internal edges" if not internal_label_used else None
			internal_label_used = True

		axes.plot([start[0], end[0]], [start[1], end[1]], color=color, linewidth=width, label=label)

	axes.set_aspect("equal", adjustable="box")
	axes.set_title("DXF external and internal edges")
	axes.set_xlabel("X")
	axes.set_ylabel("Y")
	axes.grid(False)
	axes.legend(loc="upper right")
	plt.tight_layout()
	plt.savefig(OUTPUT_IMAGE_PATH, dpi=220, bbox_inches="tight")
	plt.show()


def main() -> None:
	dxf_file_path = parse_arguments()
	document = ezdxf.readfile(str(dxf_file_path))
	# modelspace = load_modelspace(file_path=dxf_file_path)
	modelspace = document.modelspace()
	edges = collect_edges(modelspace=modelspace)

	external_count = sum(1 for edge in edges if edge[2])
	internal_count = len(edges) - external_count

	print(f"External edges: {external_count}")
	print(f"Internal edges: {internal_count}")

	plot_edges(edges=edges)


if __name__ == "__main__":
	main()
