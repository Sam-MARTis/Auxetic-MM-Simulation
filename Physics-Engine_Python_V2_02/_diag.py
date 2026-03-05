import os, sys
os.environ["SDL_VIDEODRIVER"] = "dummy"
os.environ["SDL_AUDIODRIVER"] = "dummy"

import warnings
warnings.simplefilter("always")
import numpy as np
np.set_printoptions(linewidth=120)

from Objects import Solid, Node, Edge, Vector
from constants import *
from Aux import (get_bounds_and_lines_and_hinge_points,
                 initialize_nodes_uniform_with_noise,
                 get_constraint_indices_and_constraint_functions)

min_distance, max_distance, n = 1e-6, 5e-2, 2
file_pattern = (f"/home/odin/Academic/Research/Metamaterials/Auxetic Mechanical MM/"
                f"Bistable/Computational Exploration/Modules/pattern_{n}x{n}.dxf")
bounds, lines_vec, hinge_points = get_bounds_and_lines_and_hinge_points(
    filename=file_pattern, min_distance=min_distance,
    max_distance=max_distance, center=True, size=WORLD_SIZE/2)

dxf_edges = []
for line in lines_vec:
    a, b = line[0], line[1]
    na = Node(mass=NODE_MASS, position=Vector(a.x, a.y), velocity=Vector(0,0))
    nb = Node(mass=NODE_MASS, position=Vector(b.x, b.y), velocity=Vector(0,0))
    dxf_edges.append(Edge(k=LINK_STRENGTH, l0=(b-a).magnitude(), node_a=na, node_b=nb, attach=True))

Node._count = 0
nodes, _ = initialize_nodes_uniform_with_noise(
    bounds=bounds, node_count_x=int(np.sqrt(NODE_COUNT)),
    node_count_y=int(np.sqrt(NODE_COUNT)), get_boundary_nodes=True)
hinge_nodes = [Node(mass=NODE_MASS, position=hp, velocity=Vector(0,0)) for hp in hinge_points]
nodes.extend(hinge_nodes)

solid = Solid(nodes=nodes)
solid.initialize_edges_distance_based(distance_threshold=LINK_DISTANCE)
solid.remove_intersecting_edges(dxf_edges)
solid.remove_intersecting_edges(dxf_edges)
ci, cf = get_constraint_indices_and_constraint_functions(nodes=solid.nodes, bounds=bounds)
solid.prime()
solid.set_constraint(coordinate_indices_to_constrain=ci, curve_functions=cf)

print(f"N={solid.num_nodes}, DOFs={2*solid.num_nodes}, edges={solid.num_edges}")
print(f"fixed DOFs: {len(solid.fixed_indices)}  -> {solid.fixed_indices[:10]}...")
print(f"free  DOFs: {len(solid.free_indices)}")

# ── build H exactly as solve() does ──────────────────────────────────────────
solid.cache_edges_extension()
g, H = solid.assemble_gradient_and_hessian()

print(f"\nBefore BCs:")
print(f"  rank(H) = {np.linalg.matrix_rank(H)}  (full rank = {2*solid.num_nodes})")
diag_before = np.diag(H)
print(f"  zero-diagonal DOFs: {np.sum(diag_before == 0.0)}")
print(f"  first 10 zero-diag DOFs: {np.where(diag_before == 0.0)[0][:10]}")

g, H = solid.apply_constraints_on_hessian_and_gradient(g, H)

print(f"\nAfter BCs + regularisation:")
print(f"  rank(H) = {np.linalg.matrix_rank(H)}  (full rank = {2*solid.num_nodes})")
diag_after = np.diag(H)
print(f"  zero-diagonal DOFs: {np.sum(diag_after == 0.0)}")
print(f"  min diagonal: {diag_after.min():.3e}")
print(f"  H symmetric: {np.allclose(H, H.T)}")

# ── null space analysis ───────────────────────────────────────────────────────
U, s, Vt = np.linalg.svd(H)
print(f"\nSmallest 6 singular values: {s[-6:]}")
null_mask = s < 1e-6 * s[0]
print(f"Null space dimension (tol=1e-6*smax): {null_mask.sum()}")
if null_mask.sum() > 0:
    # show which DOFs dominate each null vector
    for k in range(min(null_mask.sum(), 4)):
        v = Vt[-(k+1)]
        top = np.argsort(np.abs(v))[-6:][::-1]
        print(f"  null vec {k}: dominant DOFs {top} (node ids {top//2}), vals {v[top].round(4)}")

print("\n=== Boundary filter trace ===")
NONBOUNDED_MIDDLE_SIZE = 0.2
for node in sorted(solid.nodes, key=lambda n: n.id):
    left  = node.position.x < bounds[0]+2*EPSILON_BIG
    right = node.position.x > bounds[1]-2*EPSILON_BIG
    up    = node.position.y < bounds[2]+2*EPSILON_BIG
    down  = node.position.y > bounds[3]-2*EPSILON_BIG
    middle = np.abs(node.position.y - (bounds[2]+bounds[3])/2) < NONBOUNDED_MIDDLE_SIZE*(bounds[3]-bounds[2])
    if left or right:
        status = "SKIP" if (middle or up or down) else "FIX"
        print(f"  n{node.id:3d} x={node.position.x:7.4f} y={node.position.y:7.4f} L={int(left)} R={int(right)} up={int(up)} dn={int(down)} mid={int(middle)} -> {status}")
