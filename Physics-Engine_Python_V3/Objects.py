from typing import List, Callable
import pygame
import numpy as np
from constants import *
from multiprocessing import Pool, cpu_count
import numba as nb
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import warnings
import time
from collections import defaultdict
from scipy.spatial import cKDTree

class Helper:
    @staticmethod
    def find_lo_between_nodes(node_a: 'Node', node_b: 'Node') -> float:
        dx = node_a.position.x - node_b.position.x
        dy = node_a.position.y - node_b.position.y
        return (dx ** 2 + dy ** 2) ** 0.5
    def find_dist_sq_between_nodes(node_a: 'Node', node_b: 'Node') -> float:
        dx = node_a.position.x - node_b.position.x
        dy = node_a.position.y - node_b.position.y
        return dx ** 2 + dy ** 2


class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def dot(self, other: 'Vector') -> float:
        return self.x * other.x + self.y * other.y
    def magnitude(self) -> float:
        return (self.x ** 2 + self.y ** 2) ** 0.5
    def __sub__(self, other: 'Vector') -> 'Vector':
        return Vector(self.x - other.x, self.y - other.y)
    def __add__(self, other: 'Vector') -> 'Vector':
        return Vector(self.x + other.x, self.y + other.y)
    def __truediv__(self, scalar: float) -> 'Vector':
        a = 1/scalar
        return Vector(self.x * a, self.y * a)
    def __mul__(self, scalar: float) -> 'Vector':
        return Vector(self.x * scalar, self.y * scalar)
    def __rmul__(self, scalar: float) -> 'Vector':
        return self * scalar
    def __repr__(self):
        return f"v̄({self.x:.5g}, {self.y:.5g})"
    def make_pygame_vector(self) -> pygame.Vector2:
        return pygame.Vector2(self.x, self.y)

def world_to_screen(position: Vector, world_size:float, magnification:float) -> Vector:
    return Vector((position.x+world_size/2)*magnification , (-(position.y)+world_size/2)*magnification + PLOT_H/2)
def world_to_screen_02(position: np.ndarray, world_size:float, magnification:float) -> np.ndarray:
    x = (position[0] + world_size / 2) * magnification
    y = (-(position[1]) + world_size / 2) * magnification + PLOT_H / 2
    return np.array([x, y])
def world_to_screen_03(positions:list[tuple[float, float]], world_size:float, magnification:float) -> list[tuple[float, float]]:
    screen_positions = []
    for pos in positions:
        x = (pos[0] + world_size / 2) * magnification
        y = (-(pos[1]) + world_size / 2) * magnification + PLOT_H / 2
        screen_positions.append((x, y))
    return screen_positions

class Node:
    _count = 0
    
    def __init__(self, mass: float, position: Vector, velocity: Vector):
        self.id = Node._count
        Node._count += 1
        self.connections = []
        self.position = position
        self.velocity = velocity
        self.mass = mass
        if mass==0:
            raise ValueError("Mass cannot be zero")
        self.mass_inv = 1.0 / mass
 

def orientation(p:np.ndarray, q:np.ndarray, r):
    val = (q[1]-p[1])*(r[0]-q[0]) - (q[0]-p[0])*(r[1]-q[1])
    if abs(val) < EPSILON: return 0
    return 1 if val > 0 else -1

def segments_intersect_2d_fast(p1, p2, p3, p4):
    o1 = orientation(p1, p2, p3)
    o2 = orientation(p1, p2, p4)
    o3 = orientation(p3, p4, p1)
    o4 = orientation(p3, p4, p2)
    return o1 != o2 and o3 != o4


def nearest_points_on_line_segments(a0: Vector, a1: Vector, b0: Vector, b1: Vector):
    r = b0 - a0
    u = a1 - a0
    v = b1 - b0

    ru = r.dot(u)
    rv = r.dot(v)
    uu = u.dot(u)
    uv = u.dot(v)
    vv = v.dot(v)

    det = uu * vv - uv * uv
    s, t = 0, 0

    if det < 1e-6 * uu * vv:
        s = max(0, min(1, ru / uu))
        t = 0
    else:
        s = max(0, min(1, (ru * vv - rv * uv) / det))
        t = max(0, min(1, (ru * uv - rv * uu) / det))

    S = max(0, min(1, (t * uv + ru) / uu))
    T = max(0, min(1, (s * uv - rv) / vv))

    A:Vector = a0 + u * S
    B:Vector = b0 + v * T
    return A, B, (B - A).magnitude()

class Solid:
    def __init__(self, nodes: List[Node]):
        self.nodes = nodes
        self.nodes.sort(key=lambda node: node.id)
        self.num_nodes = len(nodes)
        
        node_positions:np.ndarray = np.zeros(DIMS*self.num_nodes, dtype=np.float64)

        for i in range(self.num_nodes):
            node_positions[2*i] = self.nodes[i].position.x
            node_positions[2*i+1] = self.nodes[i].position.y
        self.node_positions = node_positions
        
        self.t = 0
        
    


    def set_constraint(self, coordinate_indices_to_constrain: list[int], curve_functions: list[Callable[[float], float]]):
        self.coordinate_indices_to_constrain = coordinate_indices_to_constrain
        self.constraint_node_mask = np.zeros(2*self.num_nodes, dtype=bool)
        for idx in coordinate_indices_to_constrain:
            self.constraint_node_mask[idx] = True
        self.curve_functions = curve_functions
        
        self.free_indices  = np.where(~self.constraint_node_mask)[0]
        self.fixed_indices = np.where( self.constraint_node_mask)[0]


    @profile_function
    def evolve_constraint(self, t:float):
        for i in range(len(self.coordinate_indices_to_constrain)):
            idx = self.coordinate_indices_to_constrain[i]
            self.node_positions[idx] = self.curve_functions[i](t)

    @profile_function
    def cache_edges_extension(self):
        pos = self.node_positions.reshape(self.num_nodes, DIMS).T  # (DIMS, num_nodes)

        node1_ids = self.node_edges_connection[0::2]  # (num_edges,)
        node2_ids = self.node_edges_connection[1::2]

        dij = pos[:, node1_ids] - pos[:, node2_ids]  # (DIMS, num_edges)
        norms = np.linalg.norm(dij, axis=0)           # (num_edges,)

        collapsed_mask = norms <= EPSILON
        self.edges_extension[:] = np.where(collapsed_mask, 0.0, norms)

        if np.any(collapsed_mask):
            idx = np.where(collapsed_mask)[0]
            warnings.warn(
                f"[cache_edges_extension] {len(idx)} collapsed edge(s) "
                f"(d <= {EPSILON}): " +
                ", ".join(
                    f"edge {i} (n{node1_ids[i]}-n{node2_ids[i]}, d={norms[i]:.2e})"
                    for i in idx
                )
            )

    @profile_function
    def assemble_gradient_and_hessian(self):
        N   = self.num_nodes
        DOF = DIMS * N
        g   = np.zeros(DOF, dtype=np.float64)

        nec  = self.node_edges_connection         
        n1   = nec[0::2]                           
        n2   = nec[1::2]                          
        k    = self.edges_properties               
        d    = self.edges_extension                
        l0   = self.edges_extension_natural        

        pos  = self.node_positions.reshape(N, DIMS).T

        live = d != 0.0
        n1, n2, k, d, l0 = n1[live], n2[live], k[live], d[live], l0[live]

        dx = pos[0, n1] - pos[0, n2]              
        dy = pos[1, n1] - pos[1, n2]              

        inv   = 1.0 / d                           
        alpha = k * (1.0 - l0 * inv)              
        beta  = k * l0 * inv                      

        fx = alpha * dx
        fy = alpha * dy
        np.add.at(g, 2*n1,   fx);  np.add.at(g, 2*n1+1, fy)
        np.add.at(g, 2*n2,  -fx);  np.add.at(g, 2*n2+1,-fy)

        ex = dx * inv                             
        ey = dy * inv                             

        c00 = alpha + beta * ex * ex
        c01 = beta * ex * ey
        c11 = alpha + beta * ey * ey

        r1x, r1y = 2*n1,   2*n1+1
        r2x, r2y = 2*n2,   2*n2+1

        rows = np.concatenate([r1x,r1x,r1y,r1y, r2x,r2x,r2y,r2y,
                                r1x,r1x,r1y,r1y, r2x,r2x,r2y,r2y])
        cols = np.concatenate([r1x,r1y,r1x,r1y, r2x,r2y,r2x,r2y,
                                r2x,r2y,r2x,r2y, r1x,r1y,r1x,r1y])
        vals = np.concatenate([ c00, c01, c01, c11,
                                 c00, c01, c01, c11,
                                -c00,-c01,-c01,-c11,
                                -c00,-c01,-c01,-c11])

        H = sp.csc_matrix((vals, (rows, cols)), shape=(DOF, DOF))
        return g, H

    @profile_function
    def apply_constraints_on_hessian_and_gradient(self, g: np.ndarray, H: sp.csc_matrix):

        H = H.tolil() 
        fixed = self.fixed_indices
        for idx in fixed:
            H.rows[idx] = [idx]
            H.data[idx] = [1.0]
        H = H.tocsc()
    
        for idx in fixed:
            start, end = H.indptr[idx], H.indptr[idx + 1]
            H.data[start:end] = 0.0

            diag_pos = np.searchsorted(H.indices[start:end], idx)
            if diag_pos < (end - start) and H.indices[start + diag_pos] == idx:
                H.data[start + diag_pos] = 1.0
        
        g[fixed] = 0.0
        
        return g, H

    @profile_function
    def solve(self, tol, max_iter):


        for iteration in range(max_iter):
            self.cache_edges_extension()
            g, H = self.assemble_gradient_and_hessian()
            g, H = self.apply_constraints_on_hessian_and_gradient(g, H)
            dx = spsolve(H, -g)
            self.node_positions += dx*RELAXATION
            for i, idx in enumerate(self.coordinate_indices_to_constrain):
                self.node_positions[idx] = self.curve_functions[i](self.t)

            free_norm = np.linalg.norm(dx[self.free_indices])
            if free_norm < tol:
                break
        else:
            print(f"Warning: did not converge in {max_iter} iterations (||dx||={free_norm:.3e})")

    @profile_function
    def initialize_edges_distance_based(self, distance_threshold, edges_external):

        coords = np.reshape(self.node_positions, (self.num_nodes, 2))
        
        tree = cKDTree(coords)
        pairs = tree.query_pairs(distance_threshold)
        
        edges = np.zeros((len(pairs)*2), dtype=int)
        num_external_edges = len(edges_external)
        # external_edge_positions = []
        filtered_pairs = []
        # for edge in edges_external:
        #     node_a = edge.node_a
        #     node_b = edge.node_b
        #     external_edge_positions.append([np.array([node_a.position.x, node_a.position.y]), np.array([node_b.position.x, node_b.position.y])])
        for i, j in pairs:
            node_a_pos = coords[i]
            node_b_pos = coords[j]
            for k in range(num_external_edges):
                ext_node_a_pos = edges_external[k][0]
                ext_node_b_pos = edges_external[k][1]
                if segments_intersect_2d_fast(node_a_pos, node_b_pos, ext_node_a_pos, ext_node_b_pos):
                    break
            else:
                filtered_pairs.append((i, j))
  
        filtered_pairs = np.array(filtered_pairs, dtype=np.int32)

        node_1_pos = coords[filtered_pairs[:, 0]]
        node_2_pos = coords[filtered_pairs[:, 1]]
        distances = np.linalg.norm(node_1_pos - node_2_pos, axis=1)
        self.edges_properties = LINK_STRENGTH * LINK_DISTANCE / distances
        self.edges_extension_natural = distances
        self.node_edges_connection = filtered_pairs.flatten().astype(np.int32)
        self.num_edges = len(self.node_edges_connection) // 2
        self.edges_extension = np.zeros(self.num_edges, dtype=np.float64)
    
        return self.node_edges_connection



    def compute_energy_from_extension_cache(self):
        self.energy = 0.0
        self.energy = 0.5*np.sum(self.edges_properties * (self.edges_extension - self.edges_extension_natural) ** 2)
        

    def step(self, substeps:int, report_energy:bool):
        for _ in range(substeps):
            self.evolve_constraint(self.t)
            self.solve(tol=TOLERANCE_SOLVER, max_iter=MAX_ITER_SOLVER)
        if report_energy:
            self.compute_energy_from_extension_cache()
            return self.energy