from Objects import Vector, Node, Edge, Helper
from typing import Callable, List, Tuple
from constants import *
import numpy as np
import pygame



def initialize_nodes(node_count) -> list[Node]:
    nodes: list[Node] = []
    x_values = np.random.uniform(-WORLD_SIZE/4, WORLD_SIZE/4, int(node_count))
    y_values = np.random.uniform(-WORLD_SIZE/4, WORLD_SIZE/4, int(node_count))
    vx_values = np.random.uniform(-STARTING_VEL, STARTING_VEL, int(node_count))
    vy_values = np.random.uniform(-STARTING_VEL, STARTING_VEL, int(node_count))
    for x, y, vx, vy in zip(x_values, y_values, vx_values, vy_values):
        nodes.append(Node(mass=NODE_MASS, position=Vector(x, y), velocity=Vector(vx, vy)))
    return nodes

def edge_already_exists(node_a: Node, node_b: Node) -> bool:
    node_a_edges = [edge.id for edge in node_a.edges]
    node_b_edges = [edge.id for edge in node_b.edges]
    return any(edge in node_b_edges for edge in node_a_edges)

def initialize_edges_random(count, nodes: List[Node]) -> List[Edge]:
    edges: List[Edge] = []
    for i in range(count):
        node_a = nodes[np.random.randint(len(nodes))]
        node_b = nodes[np.random.randint(len(nodes))]
        if node_a == node_b:
            continue
        if edge_already_exists(node_a, node_b):
            continue
        edge = Edge(k=100.0, l0=0.8, node_a=node_a, node_b=node_b, attach=True)
        edges.append(edge)
    return edges

def initialize_edges_distance_based(count, nodes: List[Node]) -> List[Edge]:
    edges: List[Edge] = []
    num_nodes = count
    for i in range(num_nodes):
        for j in range(i+1, num_nodes):
            node_a = nodes[i]
            node_b = nodes[j]
            distance_sq = Helper.find_dist_sq_between_nodes(node_a, node_b)
            if distance_sq < LINK_DISTANCE ** 2:
                edge = Edge(k=100.0, l0=(distance_sq**0.5), node_a=node_a, node_b=node_b, attach=True)
                edges.append(edge)
    return edges

def update_edges_coordinate_properties(edges: List[Edge]):
    for edge in edges:
        edge.update_line_properties()

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

def all_edges_intersection_check(edges: List[Edge]):
    update_edges_coordinate_properties(edges)
    num_edges = len(edges)
    for i in range(num_edges):
        for j in range(i + 1, num_edges):
            e1 = edges[i]
            e2 = edges[j]
            if is_intersecting(e1, e2):
                e1.is_intersecting = True
                e2.is_intersecting = True
            # else:
            #     e1.is_intersecting = False
            #     e2.is_intersecting = False

def step_nodes(nodes: List[Node], dt: float):
	for node in nodes:
		node.step(dt)
def apply_edge_forces(edges: List[Edge], dt: float):
	for edge in edges:
		edge.step(dt)

def apply_gravity(nodes: List[Node], dt: float):
	for node in nodes:
		node.apply_force(Vector(0, -GRAVITY * node.mass), dt)




def render_everything(nodes: List[Node], edges: List[Edge], screen: pygame.Surface):
    for edge in edges:
        edge.render(screen=screen, colour=EDGE_COLOUR_SELECTED if edge.is_intersecting else EDGE_COLOUR, width=EDGE_WIDTH)
        edge.is_intersecting = False
    for node in nodes:
        node.render(screen=screen, radius=NODE_RADIUS, colour=NODE_COLOUR)