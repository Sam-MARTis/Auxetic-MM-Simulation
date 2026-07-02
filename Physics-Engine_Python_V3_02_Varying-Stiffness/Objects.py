import pygame
import warnings
import csv
import subprocess
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.sparse.linalg import spsolve
from typing import List, Callable
from Constants import *
from scipy.optimize import minimize

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

def world_to_screen_04(positions_list: np.ndarray[np.ndarray[float, float]], world_size:float, magnification:float) -> np.ndarray:
    mapped_positions = []
    for position in positions_list:
        x = (position[0] + world_size / 2) * magnification
        y = (-(position[1]) + world_size / 2) * magnification + PLOT_H / 2
        mapped_positions.append((x, y))
    return mapped_positions

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
    def __init__(self, nodes: List[Node], bounds: tuple[float, float, float, float]):
        self.nodes = nodes
        self.bounds = bounds
        self.nodes.sort(key=lambda node: node.id)
        self.num_nodes = len(nodes)
        
        node_positions:np.ndarray = np.zeros(DIMS*self.num_nodes, dtype=np.float64)

        for i in range(self.num_nodes):
            node_positions[2*i] = self.nodes[i].position.x
            node_positions[2*i+1] = self.nodes[i].position.y
        self.node_positions = node_positions
        
        self.energy = 0.0
        self.t = 0
        self.t_log = []
        self.energy_log = []
        self.charged_nodes_positions_log = []
        
    


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
    def assemble_electrostatic_gradient_and_hessian(self):
        N   = self.num_nodes
        DOF = DIMS * N
        g   = np.zeros(DOF, dtype=np.float64)

        ci  = self.charged_indices                        # (m,)
        pos = self.node_positions.reshape(N, DIMS).T      # (2, N)
        Q   = self.charges                                # (m,) charge magnitudes

        # all pairs of charged nodes
        idx_i, idx_j = np.triu_indices(len(ci), k=1)     # upper triangle, no self
        n1 = ci[idx_i]                                    # global node indices
        n2 = ci[idx_j]

        dx  = pos[0, n1] - pos[0, n2]
        dy  = pos[1, n1] - pos[1, n2]
        d   = np.sqrt(dx**2 + dy**2)

        live = d > EPSILON
        n1, n2, dx, dy, d = n1[live], n2[live], dx[live], dy[live], d[live]
        q1q2 = Q[idx_i[live]] * Q[idx_j[live]]

        inv  = 1.0 / d
        inv2 = inv * inv
        inv3 = inv2 * inv
        ex   = dx * inv
        ey   = dy * inv

        # Gradient: -Q1Q2/d^2 * e_hat
        gx = -q1q2 * inv2 * ex
        gy = -q1q2 * inv2 * ey
        np.add.at(g, 2*n1,   gx);  np.add.at(g, 2*n1+1, gy)
        np.add.at(g, 2*n2,  -gx);  np.add.at(g, 2*n2+1,-gy)

        # Hessian blocks: Q1Q2/d^3 * (I - 3*e*e^T)
        scale = q1q2 * inv3
        c00 =  scale * (1.0 - 3.0 * ex*ex)
        c01 =  scale * (     -3.0 * ex*ey)
        c11 =  scale * (1.0 - 3.0 * ey*ey)

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
    def assemble_spring_gradient_and_hessian(self):
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
    def assemble_gradient(self, node_positions: np.ndarray, edges_extension: np.ndarray, edges_properties: np.ndarray):
        N   = self.num_nodes
        DOF = DIMS * N
        g   = np.zeros(DOF, dtype=np.float64)

        nec  = self.node_edges_connection         
        n1   = nec[0::2]                           
        n2   = nec[1::2]                          
        k    = edges_properties               
        d    = edges_extension                
        l0   = self.edges_extension_natural        

        pos  = node_positions.reshape(N, DIMS).T

        live = d != 0.0
        n1, n2, k, d, l0 = n1[live], n2[live], k[live], d[live], l0[live]

        dx = pos[0, n1] - pos[0, n2]              
        dy = pos[1, n1] - pos[1, n2]              

        inv   = 1.0 / d                           
        alpha = k * (1.0 - l0 * inv)              
                          

        fx = alpha * dx
        fy = alpha * dy
        np.add.at(g, 2*n1,   fx);  np.add.at(g, 2*n1+1, fy)
        np.add.at(g, 2*n2,  -fx);  np.add.at(g, 2*n2+1,-fy)
        return g

    def assemble_gradient_and_hessian(self):
        g_spring, H_spring = self.assemble_spring_gradient_and_hessian()
        g_charge,  H_charge  = self.assemble_electrostatic_gradient_and_hessian()
        g_total = g_spring + g_charge
        H_total = H_spring + H_charge
        return g_total, H_total
    
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
    def solve(self, tol, max_iter, relaxation):


        
        for iteration in range(max_iter):
            
            self.cache_edges_extension()
            g, H = self.assemble_spring_gradient_and_hessian()
            g, H = self.apply_constraints_on_hessian_and_gradient(g, H)
            dx = spsolve(H, -g)
            self.node_positions += dx*relaxation
            # for i, idx in enumerate(self.coordinate_indices_to_constrain):
            #     self.node_positions[idx] = self.curve_functions[i](self.t)

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
        self.edges_properties = (LINK_STRENGTH * LINK_DISTANCE / distances) * np.ones(len(filtered_pairs), dtype=np.float64)
        self.edges_properties_original = self.edges_properties.copy()
        self.edges_extension_natural = distances
        self.node_edges_connection = filtered_pairs.flatten().astype(np.int32)
        self.num_edges = len(self.node_edges_connection) // 2
        self.edges_extension = np.zeros(self.num_edges, dtype=np.float64)
    
        return self.node_edges_connection



    def compute_and_cache_strain_energy_from_extension_cache(self):
        self.energy = 0.0
        self.energy = 0.5*np.sum(self.edges_properties * (self.edges_extension - self.edges_extension_natural) ** 2)
        return self.energy

    def compute_A(self, pos:np.ndarray):

        diff = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(diff, axis=2)

        with np.errstate(divide='ignore', invalid='ignore'):
            A = np.where(dist > 0, 1.0 / dist, 0.0)
        np.fill_diagonal(A, 0.0)
        return A
    
    # def add_charge_energy(self):
    #     node_positions = self.node_positions.reshape(self.num_nodes, DIMS)
    #     charged_positions = node_positions[self.charged_indices]
    #     A = self.compute_A(charged_positions)
    #     q = self.charges
    #     charge_energy = 0.5 * q @ A @ q
    #     self.energy += charge_energy
    #     print(f"Charge energy: {charge_energy:.5g}")
    #     # return charge_energy
        


    # def assign_charged_indices(self, charged_indices:np.ndarray):
    #     self.charged_indices = charged_indices
    #     self.num_charges = len(charged_indices)
    #     self.charges = np.zeros(self.num_charges, dtype=np.float64)
        

    def log_time(self):
        self.t_log.append(self.t)
    def log_energy(self):
        self.energy_log.append(self.energy)
    # def log_charged_nodes_positions(self):
    #     node_positions = self.node_positions.reshape(self.num_nodes, DIMS)
    #     self.charged_nodes_positions_log.append(node_positions[self.charged_indices].copy())

    def log_state(self):
        self.log_time()
        self.log_energy()
        # self.log_charged_nodes_positions()

    # def render_charged_points_trajectory(self, screen):
    
    #     trajectories = np.array(self.charged_nodes_positions_log).transpose(1, 0, 2)
    #     for trajectory in trajectories:
    #         screen_positions = world_to_screen_04(trajectory, world_size=WORLD_SIZE, magnification=MAGNIFICATION)
    #         pygame.draw.lines(surface=screen, color=TRACKING_COLOUR, points=screen_positions, closed=False, width=TRACKING_EDGE_WIDTH)
    
    def write_data_to_csv(self, filename):
        headings = ['time', 'strain energy'] + [f'charge_{i}_x' for i in range(self.num_charges)] + [f'charge_{i}_y' for i in range(self.num_charges)]
        assert len(self.t_log) == len(self.energy_log) == len(self.charged_nodes_positions_log), "Logs must have the same length"
        num_logs = len(self.t_log)
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(headings)
            for i in range(num_logs):
                row = [self.t_log[i], self.energy_log[i]]
                charged_positions = self.charged_nodes_positions_log[i]
                for pos in charged_positions:
                    row.extend(pos)
                writer.writerow(row)
    
    
    def step(self, substeps:int, report_energy:bool, relaxation:float):
        self.evolve_constraint(self.t)
        for _ in range(substeps):
            self.solve(tol=TOLERANCE_SOLVER, max_iter=MAX_ITER_SOLVER, relaxation=relaxation)
            self.compute_and_cache_strain_energy_from_extension_cache()

            # self.add_charge_energy()
        if report_energy:
            return self.energy
        
    def optimize_charges(self, worker_file:str, position_storage_file:str, optimized_charge_storage_file:str, verbose:bool):
        subprocess.run(["python", worker_file, position_storage_file, optimized_charge_storage_file, str(verbose)], check=True)
        with open(optimized_charge_storage_file, mode='r') as file:
            #             Charge
            # 0.12361509150398144
            # -0.4453967388332859
            # -1.2084913415454135
            # 0.35135028263123924
            # -0.6788568714239562
            # 0.9884685125697185
            # 0.733618112261491
            # -0.7055346504028998
            # 0.9031651207564099
            # 0.3452391928584101
            reader = csv.reader(file)
            next(reader)
            optimized_charges = []
            for row in reader:
                optimized_charges.append(float(row[0]))

        self.charges = np.array(optimized_charges, dtype=np.float64)* CHARGE_CONSTANT
    @profile_function
    def compute_edge_properties_from_field(self, node_positions:np.ndarray, G:np.ndarray, λx:float, λy:float, α0:float):
        inv_λx = 1.0 / λx
        inv_λy = 1.0 / λy
        num_coefficients = len(G) // 4
        
        # Extract node indices for all edges at once
        n1_ids = self.node_edges_connection[0::2]  # (num_edges,)
        n2_ids = self.node_edges_connection[1::2]
        
        # Reshape and extract node positions for all edges
        node_pos_2d = node_positions.reshape(-1, 2)  # (num_nodes, 2)
        node_a_pos = node_pos_2d[n1_ids]            # (num_edges, 2)
        node_b_pos = node_pos_2d[n2_ids]
        
        # Compute edge centroids for all edges at once
        centroids = (node_a_pos + node_b_pos) * 0.5  # (num_edges, 2)
        cx = centroids[:, 0]  # (num_edges,)
        cy = centroids[:, 1]
        max_cx = np.max(np.abs(cx))
        max_cy = np.max(np.abs(cy))
        min_cx = np.min(np.abs(cx))
        min_cy = np.min(np.abs(cy))
        
        # Precompute basis function frequencies
        freqs = np.arange(1, num_coefficients + 1, dtype=np.float64) * np.pi  # (num_coefficients,)
        
        # Compute basis functions using broadcasting: (num_edges, num_coefficients)
        sin_x = np.sin(freqs[np.newaxis, :] * cx[:, np.newaxis] * inv_λx)
        cos_x = np.cos(freqs[np.newaxis, :] * cy[:, np.newaxis] * inv_λx)
        sin_y = np.sin(freqs[np.newaxis, :] * cx[:, np.newaxis] * inv_λy)
        cos_y = np.cos(freqs[np.newaxis, :] * cy[:, np.newaxis] * inv_λy)
        
        # Reshape G into 4 coefficient groups for clarity
        G_reshaped = G.reshape(4, num_coefficients)  # (4, num_coefficients)
        
        # Compute K values via matrix multiplication: (num_edges, num_coefficients) @ (num_coefficients,) = (num_edges,)
        K1 = sin_x @ G_reshaped[0, :]
        K2 = cos_x @ G_reshaped[1, :]
        K3 = sin_y @ G_reshaped[2, :]
        K4 = cos_y @ G_reshaped[3, :]
        
        # Compute final properties for all edges
        properties = self.edges_properties_original * (1.0 + α0 * (K1 + K2 + K3 + K4))
        
        return properties
    def update_edge_properties(self, G:np.ndarray, λx:float, λy:float, α0:float):
        # inv_λx = 1.0/λx
        # inv_λy = 1.0/λy
        # num_coefficients = G.shape[1]
        # for i in range(self.num_edges):
        #     node_a_id = self.node_edges_connection[2*i]
        #     node_b_id = self.node_edges_connection[2*i+1]
        #     node_a_pos = self.node_positions[2*node_a_id:2*node_a_id+2]
        #     node_b_pos = self.node_positions[2*node_b_id:2*node_b_id+2]
        #     edge_centroid = (node_a_pos + node_b_pos) / 2
        #     K1 = np.sum(G[0][j]*np.sin((j+1)*np.pi*edge_centroid[0]*inv_λx) for j in range(num_coefficients))
        #     K2 = np.sum(G[1][j]*np.cos((j+1)*np.pi*edge_centroid[1]*inv_λx) for j in range(num_coefficients))
        #     K3 = np.sum(G[2][j]*np.sin((j+1)*np.pi*edge_centroid[0]*inv_λy) for j in range(num_coefficients))
        #     K4 = np.sum(G[3][j]*np.cos((j+1)*np.pi*edge_centroid[1]*inv_λy) for j in range(num_coefficients))
        #     self.edges_properties[i] = self.edges_properties_original[i] * (1 + α0*(K1 + K2 + K3 + K4))
        properties = self.compute_edge_properties_from_field(G=G, node_positions=self.node_positions, λx=λx, λy=λy, α0=α0)
        self.edges_properties[:] = properties[:]


    def perturb_back_check(self, perturbation_dt:float, relaxation:float, t):
        t_star = t - perturbation_dt
        current_energy = self.energy
        current_nodes_positions = self.node_positions.copy()
        current_extension = self.edges_extension.copy()
        self.evolve_constraint(self.t - perturbation_dt)
        self.solve(tol=TOLERANCE_SOLVER, max_iter=MAX_ITER_SOLVER, relaxation=relaxation)
        self.compute_and_cache_strain_energy_from_extension_cache()
        new_energy = self.energy
        new_nodes_positions = self.node_positions.copy()
        new_edges_extension = self.edges_extension.copy()   
        self.node_positions[:] = current_nodes_positions[:]
        self.edges_extension[:] = current_extension[:]
        self.energy = current_energy
        return new_energy, new_nodes_positions, new_edges_extension
    def obtain_size_at_time(self, t:float):
        return np.abs(self.bounds[0] - self.bounds[1]) + 2*PARAMETER_MULTIPLIER*t
    
    def optimize_field_rewind(self, perturbation_dt:float, relaxation:float, G:np.ndarray, t:float):
        time = self.t
        current_energy = self.energy
        current_nodes_positions = self.node_positions.copy()
        current_extension = self.edges_extension.copy()
        # K_current = self.compute_edge_properties_from_field(current_nodes_positions, G, λx=self.obtain_size_at_time(time), λy=self.obtain_size_at_time(time), α0=ALPHA_0)
    
        new_energy, new_nodes_positions, new_edges_extension = self.perturb_back_check(perturbation_dt=perturbation_dt, relaxation=relaxation, t=t)
        if new_energy < current_energy:
            current_energy = new_energy
            self.node_positions[:] = new_nodes_positions[:]
            self.edges_extension[:] = new_edges_extension[:]
            self.energy = new_energy
            print(f"Energy lowered by perturbation: {current_energy:.5g} -> {new_energy:.5g}. Stepping back")
            return True
        else:
            print(f"Perturbation did not lower energy: {current_energy:.5g} -> {new_energy:.5g}. Keeping current state")

            λx_new_point = self.obtain_size_at_time(time - perturbation_dt)
            λy_new_point = self.obtain_size_at_time(time - perturbation_dt)
            λx_old_point = self.obtain_size_at_time(time)
            λy_old_point = self.obtain_size_at_time(time)
            prev_G = G.copy()
            dx = new_nodes_positions - current_nodes_positions
            p = dx/(np.linalg.norm(dx) + 1e-8)
            @profile_function
            def objective(G, prev_G, a1, a2, a3, a4):
                K_star_new_point = self.compute_edge_properties_from_field(new_nodes_positions, G, λx=λx_new_point, λy=λy_new_point, α0=ALPHA_0)
                g_star_new_point = self.assemble_gradient(node_positions=new_nodes_positions, edges_extension=new_edges_extension, edges_properties=K_star_new_point)
                K_star_old_point = self.compute_edge_properties_from_field(current_nodes_positions, G, λx=λx_old_point, λy=λy_old_point, α0=ALPHA_0)
                g_star_old_point = self.assemble_gradient(node_positions=current_nodes_positions, edges_extension=current_extension, edges_properties=K_star_old_point)
                diff_G = G - prev_G
                obj_val = a1*np.linalg.norm(g_star_new_point)**2 - a2*np.dot(g_star_old_point, new_nodes_positions - current_nodes_positions) + a3 *np.linalg.norm(diff_G)**2 + a4*np.linalg.norm(G)**2
                return obj_val
            # def objective(G, prev_G, a1, a2, a3):    
            #     K_star_old_point = self.compute_edge_properties_from_field(current_nodes_positions, G, λx=λx_old_point, λy=λy_old_point, α0=ALPHA_0)
            #     g_star_old_point = self.assemble_gradient(node_positions=current_nodes_positions, edges_extension=current_extension, edges_properties=K_star_old_point)
            #     g_star_old_point_normalized = g_star_old_point / (np.linalg.norm(g_star_old_point) + 1e-8)

            #     k1 = g_star_old_point_normalized + p
            #     obj = a1*np.dot(k1, k1) + a2*np.linalg.norm(G)**2 + a3*np.linalg.norm(g_star_old_point)**2
                # diff_G = G - prev_G
                # dx = current_nodes_positions - new_nodes_positions
                # obj_val = 0
                return obj
            a1, a2, a3, a4 = 100, 0.01, 0.1, 10000.9
            objective_value_before = objective(G, prev_G, a1, a2, a3, a4)
            print(f"Objective value before optimization: {objective_value_before:.5g}")
            print("Preparing to optimize field corefficients")
            res  = minimize(objective, G, args=(prev_G, a1, a2, a3, a4), method='L-BFGS-B')
            print(f"Field optimization success: {res.success}, message: {res.message}, final objective value: {res.fun:.5g}")
            G = res.x

            g_star_old_point = self.assemble_gradient(node_positions=current_nodes_positions, edges_extension=current_extension, edges_properties=self.compute_edge_properties_from_field(current_nodes_positions, prev_G, λx=λx_old_point, λy=λy_old_point, α0=ALPHA_0))
            g_star_old_point_normalized = g_star_old_point / (np.linalg.norm(g_star_old_point) + 1e-8)
            k1 = g_star_old_point_normalized + p
            print(f"k1 norm: {np.linalg.norm(k1):.5g}, g_star_old_point norm: {np.linalg.norm(g_star_old_point):.5g}, G norm: {np.linalg.norm(G):.5g}")

            print(f"Updated G coefficients: {G}")
            K = self.compute_edge_properties_from_field(current_nodes_positions, G, λx=self.obtain_size_at_time(time), λy=self.obtain_size_at_time(time), α0=ALPHA_0)
            self.edges_properties[:] = K[:]
            self.cache_edges_extension()
            self.compute_and_cache_strain_energy_from_extension_cache()
            new_K_energy = self.energy
            print(f"Energy after field update: {new_K_energy:.5g} (was {current_energy:.5g} before field optimization)")
            self.evolve_constraint(t=(time - perturbation_dt))
            self.solve(tol=TOLERANCE_SOLVER, max_iter=MAX_ITER_SOLVER, relaxation=relaxation)
            self.compute_and_cache_strain_energy_from_extension_cache()
            if self.energy < new_K_energy:
                print(f"Success, energy improved : {new_K_energy:.5g} -> {self.energy:.5g}")
            else:
                print(f"Failed to improve energy after optimization : {new_K_energy:.5g} -> {self.energy:.5g}. Reverting to previous state")
            return False
    def relax_field(self, G:np.ndarray, t:float, α0:float, iterations:int, relaxation:float):
        G_initial = G.copy()
        for i in range(iterations):
            factor = 1.0 - (i+1)/iterations
            G = G_initial * factor
            K = self.compute_edge_properties_from_field(self.node_positions, G, λx=self.obtain_size_at_time(t), λy=self.obtain_size_at_time(t), α0=α0)
            self.edges_properties[:] = K[:] 
            self.evolve_constraint(t=t)
            self.solve(tol=TOLERANCE_SOLVER, max_iter=MAX_ITER_SOLVER, relaxation=relaxation)
            self.compute_and_cache_strain_energy_from_extension_cache()
            





