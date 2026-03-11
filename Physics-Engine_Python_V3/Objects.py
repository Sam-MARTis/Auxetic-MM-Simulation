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


class Node:
    _count = 0
    
    def __init__(self, mass: float, position: Vector, velocity: Vector):
        self.id = Node._count
        Node._count += 1
        self.connections = []
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.edges: list[Edge] = []
        if mass==0:
            raise ValueError("Mass cannot be zero")
        self.mass_inv = 1.0 / mass
    def render(self, screen: pygame.Surface, radius: float, colour: tuple):
        screen_position = world_to_screen(self.position, WORLD_SIZE, MAGNIFICATION)
        pygame.draw.circle(screen, colour, (int(screen_position.x), int(screen_position.y)), radius)
    def step(self, dt: float):
        self.position.x += self.velocity.x * dt
        self.position.y += self.velocity.y * dt
        if abs(self.position.x) > BOUNDARY_X_P:
            self.position.x = max(min(self.position.x, BOUNDARY_X_P), BOUNDARY_X_N)
            self.velocity.x = 0
        if abs(self.position.y) > BOUNDARY_Y_P:
            self.position.y = max(min(self.position.y, BOUNDARY_Y_P), BOUNDARY_Y_N)
            self.velocity.y = 0
    def apply_force(self, force:Vector, dt:float):
        impulse =  force * (dt * self.mass_inv)
        self.velocity.x += impulse.x
        self.velocity.y += impulse.y
        self.velocity.x -= self.velocity.x * dt*VISCOCITY
        self.velocity.y -= self.velocity.y * dt*VISCOCITY
    def attach_edge(self, edge: 'Edge'):
        self.edges.append(edge)
        return self

class Edge:
    _count = 0
    def __init__(self, k: float, l0: float, node_a: Node, node_b: Node, attach:bool):
        self.k = k
        self.l0 = l0
        self.node_a = node_a
        self.node_b = node_b
        self.id = Edge._count
        Edge._count += 1
        if attach:
            self.node_a.attach_edge(self)
            self.node_b.attach_edge(self)
        self.m = 0
        self.c = 0
        self.is_intersecting = False
        self.active = attach
        self.l0_inv = max(1.0 / l0 , 50)

    def activate_edge(self, nodea, nodeb):
        self.node_a.attach_edge(self)
        self.node_b.attach_edge(self)
        self.active = True
    def update_line_properties(self):
        dx = self.node_a.position.x - self.node_b.position.x
        self.m = (self.node_a.position.y - self.node_b.position.y) / dx if dx != 0 else float('inf')
        self.c = self.node_a.position.y - self.m * self.node_a.position.x
    def update_l0(self, l0):
        self.l0 = l0
        self.l0_inv = max(1.0 / l0 , 50)
    def step(self, dt:float):
        if self.active:
            r_vec = (self.node_a.position - self.node_b.position)
            l = r_vec.magnitude()
            dl = (l - self.l0)
            
            force_magnitude = -self.k * dl * np.abs(dl)
            force_direction = r_vec / l if l != 0 else ValueError("Zero length vector in edge step")
            self.node_a.apply_force(force=force_direction * force_magnitude, dt=dt)
            self.node_b.apply_force(force= force_direction * force_magnitude * -1, dt=dt)
        else:
            return
    def render(self, screen: pygame.Surface, colour: tuple, width: int):
        if self.active:
            pygame.draw.line(surface=screen, color=colour, start_pos=world_to_screen(self.node_a.position, WORLD_SIZE, MAGNIFICATION).make_pygame_vector(), end_pos=world_to_screen(self.node_b.position, WORLD_SIZE, MAGNIFICATION).make_pygame_vector(), width=width)
    def delete_edge(self):
        self.node_a.edges.remove(self)
        self.node_b.edges.remove(self)
        self.node_a = None
        self.node_b = None
        self.active = False
    def compute_forces(self, pos_a:Vector, pos_b:Vector):
        r_vec = pos_a - pos_b
        l = r_vec.magnitude()
        dl = (l - self.l0)

        force_magnitude = -self.k * dl
        force_direction = r_vec / l if l != 0 else ValueError("Zero length vector in edge step")
        return force_direction * force_magnitude


def is_intersecting(edge_a: Edge, edge_b: Edge) -> bool:
    dm = edge_a.m - edge_b.m
    x_meet = (edge_b.c - edge_a.c) / dm if dm != 0 else float('inf')
    y_meet = edge_a.m * x_meet + edge_a.c
    hit = True
    hit_a = x_meet <( max(edge_a.node_a.position.x, edge_a.node_b.position.x)-EPSILON) and (x_meet > min(edge_a.node_a.position.x, edge_a.node_b.position.x)+EPSILON)
    hit_b = x_meet <( max(edge_b.node_a.position.x, edge_b.node_b.position.x)-EPSILON) and (x_meet > min(edge_b.node_a.position.x, edge_b.node_b.position.x)+EPSILON)
    hit = hit_a and hit_b
    # print(hit)
    return hit


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
        
        self.edges: list[Edge] = []
        self.t = 0
    
    @profile_function
    def prime(self):

        self.node_edges_connection = np.ones(2*self.num_edges, dtype=np.int32)* (-1)
        self.edges_properties = np.zeros(self.num_edges, dtype=np.float64)
        self.edges_extension_natural = np.zeros(self.num_edges, dtype=np.float64)
        self.edges_extension = np.zeros(self.num_edges, dtype=np.float64)
        for i in range(self.num_edges):
            edge = self.edges[i]
            self.node_edges_connection[2*i] = edge.node_a.id
            self.node_edges_connection[2*i+1] = edge.node_b.id
            self.edges_properties[i] = edge.k
            self.edges_extension_natural[i] = edge.l0
            self.edges_extension[i] = 0
  

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
            self.node_positions += dx
            for i, idx in enumerate(self.coordinate_indices_to_constrain):
                self.node_positions[idx] = self.curve_functions[i](self.t)

            free_norm = np.linalg.norm(dx[self.free_indices])
            if free_norm < tol:
                break
        else:
            print(f"Warning: did not converge in {max_iter} iterations (||dx||={free_norm:.3e})")


    def add_edge(self, edge: Edge):
        self.edges.append(edge)
        
    def edge_already_exists(self, node_a: Node, node_b: Node) -> bool:
        node_a_edges = [edge.id for edge in node_a.edges]
        node_b_edges = [edge.id for edge in node_b.edges]
        return any(edge in node_b_edges for edge in node_a_edges)
    
    def edge_already_exists_id(self, id_a: int, id_b: int) -> bool:
        return self.edge_already_exists(self.nodes[id_a], self.nodes[id_b])
    
    @profile_function
    def initialize_edges_distance_based(self, distance_threshold):

        coords = np.array([[n.position.x, n.position.y] for n in self.nodes])
        tree = cKDTree(coords)
        pairs = tree.query_pairs(distance_threshold)
        for i, j in pairs:
            node_a = self.nodes[i]
            node_b = self.nodes[j]

            distance = np.linalg.norm(coords[i] - coords[j])

            edge = Edge(
                k=LINK_STRENGTH * LINK_DISTANCE / distance,
                l0=distance,
                node_a=node_a,
                node_b=node_b,
                attach=True
            )
            self.add_edge(edge)

        self.num_edges = len(self.edges)
        return self.edges

    @profile_function
    def remove_intersecting_edges(self, edges_external: list[Edge]):
        if not self.edges or not edges_external:
            return

        # Shape: (N, 3) and (M, 3)
        a0 = np.array([e.node_a.position for e in self.edges])
        a1 = np.array([e.node_b.position for e in self.edges])
        b0 = np.array([e.node_a.position for e in edges_external])
        b1 = np.array([e.node_b.position for e in edges_external])

        # Broadcast: (N, 1, 3) vs (1, M, 3) → (N, M, 3)
        d1 = a1[:, None] - a0[:, None]        # (N, 1, 3)
        d2 = b1[None] - b0[None]              # (1, M, 3)
        r  = a0[:, None] - b0[None]           # (N, M, 3)

        d11 = np.einsum('ijk,ijk->ij', d1, d1)  # (N, M)
        d22 = np.einsum('ijk,ijk->ij', d2, d2)
        d12 = np.einsum('ijk,ijk->ij', d1, d2)
        rd1 = np.einsum('ijk,ijk->ij', r,  d1)
        rd2 = np.einsum('ijk,ijk->ij', r,  d2)

        denom = d11 * d22 - d12 * d12
        parallel = np.abs(denom) < EPSILON

        t = np.where(parallel, 0.0, (rd1 * d22 - rd2 * d12) / np.where(parallel, 1.0, denom))
        u = np.where(parallel, 0.0, (rd1 * d12 - rd2 * d11) / np.where(parallel, 1.0, denom))

        t = np.clip(t, 0.0, 1.0)
        u = np.clip(u, 0.0, 1.0)

        closest = (a0[:, None] + t[..., None] * d1) - (b0[None] + u[..., None] * d2)
        dist_sq = np.einsum('ijk,ijk->ij', closest, closest)  # (N, M)

        remove_mask = np.any(dist_sq < EPSILON ** 2, axis=1)  # (N,)

        for i in reversed(range(len(self.edges))):
            if remove_mask[i]:
                self.delete_edge(self.edges[i])
    def initialize_edges_distance_based_02(self, distance_threshold:float) -> List[Edge]:
        pass

    def delete_edge(self, edge: Edge):
        edge.delete_edge()
        self.edges.remove(edge)
        self.num_edges -= 1

    def step_nodes(self, dt: float):
        for node in self.nodes:
            node.step(dt)
    def apply_edge_forces(self, dt: float):
        for edge in self.edges:
            edge.step(dt)

    def apply_gravity(self, dt: float, gravity: float):
        for node in self.nodes:
            node.apply_force(Vector(0, -gravity * node.mass), dt)
    def apply_gravity_low_level(self, dt: float, gravity: float):
        for i in range(self.num_nodes):
            self.node_velocities[2*i+1] -= gravity * dt
    def compute_energy_from_extension_cache(self):
        self.energy = 0.0
        for i in range(self.num_edges):
            k  = self.edges_properties[i]
            d  = self.edges_extension[i]   
            l0 = self.edges_extension_natural[i]
            self.energy += 0.5 * k * (d - l0) ** 2
        

    def step(self, dt:float, substeps:int, report_energy:bool):
        for _ in range(substeps):

            self.evolve_constraint(self.t)
            self.solve(tol=TOLERANCE_SOLVER, max_iter=MAX_ITER_SOLVER)

        if report_energy:
            self.compute_energy_from_extension_cache()
            return self.energy