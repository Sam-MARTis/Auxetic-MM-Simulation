from Objects import Vector, Node, Edge
from typing import Callable, List
from constants import *
import numpy as np


def initialize_nodes(node_count) -> list[Node]:
    nodes: list[Node] = []
    x_values = np.random.uniform(-WORLD_SIZE/2, WORLD_SIZE/2, int(node_count))
    y_values = np.random.uniform(-WORLD_SIZE/2, WORLD_SIZE/2, int(node_count))
    vx_values = np.random.uniform(-STARTING_VEL, STARTING_VEL, int(node_count))
    vy_values = np.random.uniform(-STARTING_VEL, STARTING_VEL, int(node_count))
    for x, y, vx, vy in zip(x_values, y_values, vx_values, vy_values):
        nodes.append(Node(mass=NODE_MASS, position=Vector(x, y), velocity=Vector(vx, vy)))
    return nodes

def edge_already_exists(node_a: Node, node_b: Node) -> bool:
    node_a_edges = [edge.id for edge in node_a.edges]
    node_b_edges = [edge.id for edge in node_b.edges]
    return any(edge in node_b_edges for edge in node_a_edges)

def initialize_edges(count, nodes: List[Node]) -> List[Edge]:
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

def step_nodes(nodes: List[Node], dt: float):
	for node in nodes:
		node.step(dt)
def apply_edge_forces(edges: List[Edge], dt: float):
	for edge in edges:
		edge.step(dt)