from Objects import Vector, Node
from constants import *
import numpy as np


def initialize_nodes(node_count) -> list[Node]:
    nodes: list[Node] = []
    x_values = np.linspace(-WORLD_SIZE/2, WORLD_SIZE/2, int(np.sqrt(node_count)))
    y_values = np.linspace(-WORLD_SIZE/2, WORLD_SIZE/2, int(np.sqrt(node_count)))
    for x in x_values:
        for y in y_values:
            nodes.append(Node(Vector(x, y)))
    return nodes