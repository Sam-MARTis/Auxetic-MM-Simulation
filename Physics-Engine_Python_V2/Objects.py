from typing import List, Callable
import pygame
import numpy as np
from constants import *
from multiprocessing import Pool, cpu_count
import numba as nb
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import warnings

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
    return Vector((position.x+world_size/2)*magnification, (-(position.y)+world_size/2)*magnification)


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
        self.edges: list[Edge] = []
        self.t = 0
    
    def prime(self):
        node_positions:np.ndarray = np.zeros(DIMS*self.num_nodes, dtype=np.float64)

        for i in range(self.num_nodes):
            node_positions[2*i] = self.nodes[i].position.x
            node_positions[2*i+1] = self.nodes[i].position.y
        self.node_positions = node_positions

        self.node_edges_connection = np.ones(2*self.num_edges, dtype=np.int32)* (-1)
        self.edges_properties = np.zeros(self.num_edges, dtype=np.float64)
        self.edges_extension_natural = np.zeros(self.num_edges, dtype=np.float64)
        self.edges_extension_inv = np.zeros(self.num_edges, dtype=np.float64)
        self.K_arr = sp.dok_array((2*self.num_nodes, 2*self.num_nodes), dtype=np.float64)

        for i in range(self.num_edges):
            edge = self.edges[i]
            self.node_edges_connection[2*i] = edge.node_a.id
            self.node_edges_connection[2*i+1] = edge.node_b.id
            self.K_arr[2*edge.node_a.id, 2*edge.node_b.id] = edge.k
            self.K_arr[2*edge.node_a.id+1, 2*edge.node_b.id+1] = edge.k
            self.K_arr[2*edge.node_b.id, 2*edge.node_a.id] = edge.k
            self.K_arr[2*edge.node_b.id+1, 2*edge.node_a.id+1] = edge.k
            self.edges_properties[i] = edge.k
            self.edges_extension_natural[i] = edge.l0
            self.edges_extension_inv[i] = 1.0 / edge.l0

        self.L = np.zeros((2*self.num_nodes, 2*self.num_nodes), dtype=np.float64)
        print(f"Performing verification check of K_ij")
        for i in range(2*self.num_nodes):
            if self.K_arr[i, i] != 0:
                raise RuntimeError(f"K[{i}, {i}] = {self.K_arr[i, i]}")
        print(f"Finished verification check of K_ij")
        
        # for i in range(2*self.num_nodes):
        #     if i % 100 == 0:
        #         print(f"Processing row {i} out of {2*self.num_nodes}")
        #     for j in range(2*self.num_nodes):
        #         self.L[i, j] = -self.K_arr[i, j]
        #         self.L[j, i] = -self.K_arr[i, j]
        #         self.L[i, i] += self.K_arr[i, j]
        #         self.L[j, j] += self.K_arr[i, j]
        self.L -= self.K_arr
        for i in range(2*self.num_nodes):
            self.L[i, i] = self.L[i, :].sum()

    

    def set_constraint(self, coordinate_indices_to_constrain: list[int], curve_functions: list[Callable[[float], float]]):
        self.coordinate_indices_to_constrain = coordinate_indices_to_constrain
        self.constraint_node_mask = np.zeros(2*self.num_nodes, dtype=bool)
        for idx in coordinate_indices_to_constrain:
            self.constraint_node_mask[idx] = True
        self.curve_functions = curve_functions
        
        self.free_indices  = np.where(~self.constraint_node_mask)[0]
        self.fixed_indices = np.where( self.constraint_node_mask)[0]

        self.L_ff = self.L[self.free_indices, :][:, self.free_indices]
        self.L_fc = self.L[self.free_indices, :][:, self.fixed_indices]

        
        
    
    def evolve_constraint(self, t:float):
        # print(f"Length of coordinate_indices_to_constrain: {len(self.coordinate_indices_to_constrain)}")
        # print(f"Length of curve_functions: {len(self.curve_functions)}")
        for i in range(len(self.coordinate_indices_to_constrain)):
            idx = self.coordinate_indices_to_constrain[i]
            self.node_positions[idx] = self.curve_functions[i](t)
          
    def cache_edges_extension_inv(self):
        pos = np.reshape(self.node_positions, (DIMS, self.num_nodes))
        collapsed = []
        for i in range(self.num_edges):
            node1_id = self.node_edges_connection[2*i]
            node2_id = self.node_edges_connection[2*i+1]
            p1 = pos[:, node1_id]
            p2 = pos[:, node2_id]
            dij = p1 - p2
            norm_dij = np.linalg.norm(dij)
            if norm_dij <= EPSILON:
                self.edges_extension_inv[i] = 0.0
                collapsed.append((i, node1_id, node2_id, norm_dij))
            else:
                self.edges_extension_inv[i] = 1.0 / norm_dij
        if collapsed:
            warnings.warn(
                f"[cache_edges_extension_inv] {len(collapsed)} collapsed edge(s) "
                f"(d <= {EPSILON}): " +
                ", ".join(f"edge {ei} (n{a}-n{b}, d={d:.2e})" for ei, a, b, d in collapsed)
            )


    def assemble_gradient_and_hessian(self):
        g = np.zeros((DIMS*self.num_nodes), dtype=np.float64)
        C = np.zeros((DIMS*self.num_nodes, DIMS*self.num_nodes), dtype=np.float64)
        pos = np.reshape(self.node_positions, (DIMS, self.num_nodes))

        for i in range(self.num_edges):
            node1_id = self.node_edges_connection[2*i]
            node2_id = self.node_edges_connection[2*i+1]
            dij = pos[:, node1_id] - pos[:, node2_id]

            k   = self.K_arr[node1_id, node2_id]
            inv = self.edges_extension_inv[i]
            l0  = self.edges_extension_natural[i]


            if inv == 0.0:
                continue

            e   = dij * inv

            f   = k * (1 - l0*inv) * dij
            g[2*node1_id:2*node1_id+2] += f
            g[2*node2_id:2*node2_id+2] -= f

            Cij = k * ((1 - l0*inv) * np.eye(DIMS) + l0*inv * np.outer(e, e))
            C[2*node1_id: 2*node1_id+2, 2*node1_id:2*node1_id+2] += Cij
            C[2*node2_id:2*node2_id+2, 2*node2_id:2*node2_id+2] += Cij
            C[2*node1_id:2*node1_id+2, 2*node2_id:2*node2_id+2] -= Cij
            C[2*node2_id:2*node2_id+2, 2*node1_id:2*node1_id+2] -= Cij
        # H = C.transpose(0, 2, 1, 3).reshape(DIMS*self.num_nodes, DIMS*self.num_nodes)
        return g, C

    def apply_constraints_on_hessian_and_gradient(self, g:np.ndarray, H:np.ndarray):
        H[self.fixed_indices, :] = 0.0
        H[:, self.fixed_indices] = 0.0
        H[self.fixed_indices, self.fixed_indices] = 1.0
        g[self.fixed_indices] = 0.0
        return g, H

    def _compute_energy(self, positions: np.ndarray) -> float:
        """E = 0.5 * sum_edges K_ij * (||x_i - x_j|| - l0)^2"""
        pos = positions.reshape(DIMS, self.num_nodes)
        E = 0.0
        for i in range(self.num_edges):
            n1 = self.node_edges_connection[2*i]
            n2 = self.node_edges_connection[2*i+1]
            dij = pos[:, n1] - pos[:, n2]
            d   = np.linalg.norm(dij)
            l0  = self.edges_extension_natural[i]
            k   = self.K_arr[n1, n2]
            E  += 0.5 * k * (d - l0) ** 2
        return E

    def solve(self, tol, max_iter):


        for iteration in range(max_iter):
            self.cache_edges_extension_inv()
            g, H = self.assemble_gradient_and_hessian()
            g, H = self.apply_constraints_on_hessian_and_gradient(g, H)

 
            print(f"Check if H is symmetric: {np.allclose(H, H.T)}")
            print(f"Check if any diagonals are zero: {np.any(np.isclose(np.diag(H), 0))}\n\n")
            H_inv = np.linalg.inv(H)
            
            dx = -H_inv @ g
    

            self.node_positions += dx*0.8

        else:
            print(f"Warning: did not converge in {max_iter} iterations")

    def add_edge(self, edge: Edge):
        self.edges.append(edge)
        
    def edge_already_exists(self, node_a: Node, node_b: Node) -> bool:
        node_a_edges = [edge.id for edge in node_a.edges]
        node_b_edges = [edge.id for edge in node_b.edges]
        return any(edge in node_b_edges for edge in node_a_edges)
    
    def edge_already_exists_id(self, id_a: int, id_b: int) -> bool:
        return self.edge_already_exists(self.nodes[id_a], self.nodes[id_b])
    
    def initialize_edges_distance_based(self, distance_threshold:float) -> List[Edge]:
        for i in range(self.num_nodes):
            for j in range(i+1, self.num_nodes):
                node_a = self.nodes[i]
                node_b = self.nodes[j]
                distance_sq = Helper.find_dist_sq_between_nodes(node_a, node_b)
                if distance_sq < distance_threshold ** 2:
                    edge = Edge(k=LINK_STRENGTH, l0=(distance_sq**0.5)*1.0, node_a=node_a, node_b=node_b, attach=True)
                    self.add_edge(edge)
        self.num_edges = len(self.edges)
        return self.edges

    def remove_intersecting_edges(self, edges_external: list[Edge]):
        remove_indices = [False] * len(self.edges)
        for i, edge in enumerate(self.edges):
            for external_edge in edges_external:
                _, _, dist = nearest_points_on_line_segments(edge.node_a.position, edge.node_b.position, external_edge.node_a.position, external_edge.node_b.position)
                if dist < EPSILON:
                    # self.delete_edge(edge)
                    remove_indices[i] = True
                    break

        for i in reversed(range(len(remove_indices))):
            if remove_indices[i]:
                self.delete_edge(self.edges[i])

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

    
    def step(self, dt:float, substeps:int):
        for _ in range(substeps):
            self.evolve_constraint(self.t)
            self.solve(tol=TOLERANCE_SOLVER, max_iter=MAX_ITER_SOLVER)
            
            self.t += dt / substeps
        pass
