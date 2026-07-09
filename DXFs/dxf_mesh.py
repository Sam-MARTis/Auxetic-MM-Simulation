"""
dxf_mesh.py
-----------
Reads a DXF file, meshes the enclosed region using Gmsh (OCC kernel),
and outputs the triangulated mesh as an image.

Usage:
    python dxf_mesh.py <input.dxf> [mesh_size] [output.png]

Arguments:
    input.dxf   : path to DXF file
    mesh_size   : target element size (default: auto, ~1/30 of bounding box diagonal)
    output.png  : output image path (default: mesh.png)
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import gmsh


def extract_boundary_from_dxf(dxf_path: str):
    """
    Extract closed polyline / spline / circle / arc boundaries from a DXF.
    Returns a list of loops, each loop being an (N,2) numpy array of points.
    """
    import ezdxf
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()

    loops = []

    for entity in msp:
        etype = entity.dxftype()

        if etype == "LWPOLYLINE":
            pts = np.array(list(entity.vertices()), dtype=float)[:, :2]
            if entity.closed or np.allclose(pts[0], pts[-1]):
                if not np.allclose(pts[0], pts[-1]):
                    pts = np.vstack([pts, pts[0]])
            loops.append(pts)

        elif etype == "POLYLINE":
            pts = np.array([v.dxf.location[:2] for v in entity.vertices], dtype=float)
            if entity.is_closed or np.allclose(pts[0], pts[-1]):
                if not np.allclose(pts[0], pts[-1]):
                    pts = np.vstack([pts, pts[0]])
            loops.append(pts)

        elif etype == "CIRCLE":
            cx, cy = entity.dxf.center.x, entity.dxf.center.y
            r = entity.dxf.radius
            t = np.linspace(0, 2 * np.pi, 64, endpoint=False)
            pts = np.column_stack([cx + r * np.cos(t), cy + r * np.sin(t)])
            pts = np.vstack([pts, pts[0]])
            loops.append(pts)

        elif etype in ("SPLINE", "ARC"):
            # Approximate via flattened points
            try:
                pts = np.array([p[:2] for p in entity.flattening(0.01)], dtype=float)
                if len(pts) >= 3:
                    loops.append(pts)
            except Exception:
                pass

    if not loops:
        raise ValueError(f"No closed boundary entities found in {dxf_path}. "
                         "Supported: LWPOLYLINE, POLYLINE, CIRCLE, SPLINE, ARC.")
    return loops


def mesh_dxf(dxf_path: str, mesh_size: float = None, output_image: str = "mesh.png"):
    # ------------------------------------------------------------------ #
    # 1. Extract boundary loops from DXF
    # ------------------------------------------------------------------ #
    loops = extract_boundary_from_dxf(dxf_path)
    print(f"Found {len(loops)} boundary loop(s).")

    # ------------------------------------------------------------------ #
    # 2. Auto mesh size from bounding box if not provided
    # ------------------------------------------------------------------ #
    all_pts = np.vstack(loops)
    xmin, ymin = all_pts.min(axis=0)
    xmax, ymax = all_pts.max(axis=0)
    diag = np.hypot(xmax - xmin, ymax - ymin)
    if mesh_size is None:
        mesh_size = diag / 30
    print(f"Bounding box: ({xmin:.3f}, {ymin:.3f}) -> ({xmax:.3f}, {ymax:.3f})")
    print(f"Mesh size   : {mesh_size:.4f}")

    # ------------------------------------------------------------------ #
    # 3. Build geometry in Gmsh (OCC kernel)
    # ------------------------------------------------------------------ #
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    gmsh.model.add("dxf_mesh")

    curve_loops = []

    for loop_idx, loop_pts in enumerate(loops):
        # Ensure loop is closed
        if not np.allclose(loop_pts[0], loop_pts[-1]):
            loop_pts = np.vstack([loop_pts, loop_pts[0]])

        # Add points
        pt_tags = []
        for x, y in loop_pts[:-1]:           # exclude repeated closing point
            tag = gmsh.model.occ.addPoint(x, y, 0.0, mesh_size)
            pt_tags.append(tag)

        # Add lines
        line_tags = []
        n = len(pt_tags)
        for i in range(n):
            tag = gmsh.model.occ.addLine(pt_tags[i], pt_tags[(i + 1) % n])
            line_tags.append(tag)

        cl = gmsh.model.occ.addCurveLoop(line_tags)
        curve_loops.append(cl)

    # First loop = outer boundary; remaining = holes
    surf = gmsh.model.occ.addPlaneSurface(curve_loops)

    gmsh.model.occ.synchronize()

    # ------------------------------------------------------------------ #
    # 4. Mesh
    # ------------------------------------------------------------------ #
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Laplace2D")

    # ------------------------------------------------------------------ #
    # 5. Extract nodes and triangles
    # ------------------------------------------------------------------ #
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = coords.reshape(-1, 3)
    x_nodes = coords[:, 0]
    y_nodes = coords[:, 1]

    # Map gmsh 1-based node tags to 0-based indices
    tag_to_idx = {tag: i for i, tag in enumerate(node_tags)}

    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(dim=2)
    triangles = []
    for et, ent in zip(elem_types, elem_node_tags):
        if et == 2:   # triangle
            tri_arr = ent.reshape(-1, 3)
            for tri in tri_arr:
                triangles.append([tag_to_idx[t] for t in tri])

    triangles = np.array(triangles)
    gmsh.finalize()

    print(f"Nodes    : {len(x_nodes)}")
    print(f"Triangles: {len(triangles)}")

    # ------------------------------------------------------------------ #
    # 6. Plot
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f"Mesh — {os.path.basename(dxf_path)}", fontsize=13)

    triang = mtri.Triangulation(x_nodes, y_nodes, triangles)

    # Left: mesh only
    ax = axes[0]
    ax.set_title(f"Triangulation  ({len(triangles)} elements)")
    ax.triplot(triang, color="steelblue", linewidth=0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Right: filled + edges, colored by element area
    ax = axes[1]
    ax.set_title("Element area distribution")
    areas = []
    for tri in triangles:
        x0, y0 = x_nodes[tri[0]], y_nodes[tri[0]]
        x1, y1 = x_nodes[tri[1]], y_nodes[tri[1]]
        x2, y2 = x_nodes[tri[2]], y_nodes[tri[2]]
        area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        areas.append(area)
    areas = np.array(areas)

    tpc = ax.tripcolor(triang, facecolors=areas, cmap="viridis", edgecolors="k",
                       linewidth=0.3)
    fig.colorbar(tpc, ax=ax, label="Element area")
    ax.set_aspect("equal")
    ax.set_xlabel("x")

    plt.tight_layout()
    plt.savefig(output_image, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_image}")
    return x_nodes, y_nodes, triangles


if __name__ == "__main__":
    dxf_path    = sys.argv[1] if len(sys.argv) > 1 else "test.dxf"
    mesh_size   = float(sys.argv[2]) if len(sys.argv) > 2 else None
    output_img  = sys.argv[3] if len(sys.argv) > 3 else "mesh.png"

    mesh_dxf(dxf_path, mesh_size, output_img)
