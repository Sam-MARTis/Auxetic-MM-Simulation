from typing import List, Callable
import pygame
from constants import *



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
    def __truediv__(self, scalar: float) -> 'Vector':
        a = 1/scalar
        return Vector(self.x * a, self.y * a)
    def __mul__(self, scalar: float) -> 'Vector':
        return Vector(self.x * scalar, self.y * scalar)
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
    def step(self, dt:float):
        if self.active:
            r_vec = (self.node_a.position - self.node_b.position)
            l = r_vec.magnitude()

            force_magnitude = -self.k * (l - self.l0)
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

class Solid:
    def __init__(self, nodes: List[Node]):
        self.nodes = nodes
        self.num_nodes = len(nodes)
        self.edges: list[Edge] = []
        
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
                    edge = Edge(k=LINK_STRENGTH, l0=(distance_sq**0.5)*1.1, node_a=node_a, node_b=node_b, attach=True)
                    self.add_edge(edge)
        return self.edges

    def delete_edge(self, edge: Edge):
        edge.delete_edge()
        self.edges.remove(edge)
        
        
    def step_nodes(self, dt: float):
        for node in self.nodes:
            node.step(dt)
    def apply_edge_forces(self, dt: float):
        for edge in self.edges:
            edge.step(dt)

    def apply_gravity(self, dt: float, gravity: float):
        for node in self.nodes:
            node.apply_force(Vector(0, -gravity * node.mass), dt)

    def step(self, dt:float, substeps:int, gravity:float):
        substeps_dt = dt / substeps
        for _ in range(substeps):
            self.step_nodes(substeps_dt)
            self.apply_gravity(dt=substeps_dt, gravity=gravity)
        self.apply_edge_forces(dt=substeps_dt)


