from Objects import Vector, Node, Edge, Helper, Solid, world_to_screen
from typing import Callable, List, Tuple
from constants import *
import numpy as np
import pygame
import ezdxf as dxf

np.random.seed(42)

def initialize_nodes(node_count) -> list[Node]:
    nodes: list[Node] = []
    x_values = np.random.uniform(-WORLD_SIZE/DENSITY_FACTOR, WORLD_SIZE/DENSITY_FACTOR, int(node_count))
    y_values = np.random.uniform(-WORLD_SIZE/DENSITY_FACTOR, WORLD_SIZE/DENSITY_FACTOR, int(node_count))
    vx_values = np.random.uniform(-STARTING_VEL, STARTING_VEL, int(node_count))
    vy_values = np.random.uniform(-STARTING_VEL, STARTING_VEL, int(node_count))
    for x, y, vx, vy in zip(x_values, y_values, vx_values, vy_values):
        nodes.append(Node(mass=NODE_MASS, position=Vector(x, y), velocity=Vector(vx, vy)))
    return nodes

def initialize_nodes_uniform(bounds: Tuple[float, float, float, float], node_count_x: int, node_count_y: int) -> list[Node]:
    nodes: list[Node] = []
    x_min, x_max, y_min, y_max = bounds
    x_values = np.linspace(x_min+EPSILON_BIG, x_max - EPSILON_BIG, node_count_x)
    y_values = np.linspace(y_min+EPSILON_BIG, y_max - EPSILON_BIG, node_count_y)
    for x in x_values:
        for y in y_values:
            nodes.append(Node(mass=NODE_MASS, position=Vector(x, y), velocity=Vector(0, 0)))
    return nodes
def initialize_nodes_uniform_with_noise(bounds: Tuple[float, float, float, float], node_count_x: int, node_count_y: int, get_boundary_nodes:bool = False) -> Tuple[list[Node], list[Node]]:
    nodes: list[Node] = []
    boundary_nodes:list[Node] = []
    x_min, x_max, y_min, y_max = bounds
    x_values = np.linspace(x_min+EPSILON_BIG, x_max - EPSILON_BIG, node_count_x)
    y_values = np.linspace(y_min+EPSILON_BIG, y_max - EPSILON_BIG, node_count_y)
    for i in range(1, len(x_values)-1):
        for j in range(1, len(y_values)-1):
            x_values[i] += np.random.normal(0, NOISE)
            y_values[j] += np.random.normal(0, NOISE)

    for x in x_values:
        for y in y_values:
            node = Node(mass=NODE_MASS, position=Vector(x, y), velocity=Vector(0, 0))
            nodes.append(node)
            if(x <= x_min + EPSILON_BIG or x >= x_max - EPSILON_BIG or y <= y_min + EPSILON_BIG or y >= y_max - EPSILON_BIG):
                boundary_nodes.append(node)
    return nodes, boundary_nodes
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
                edge = Edge(k=LINK_STRENGTH, l0=(distance_sq**0.5)*1.1, node_a=node_a, node_b=node_b, attach=True)
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
        e1 = edges[i]
        for j in range(i + 1, num_edges):
            e2 = edges[j]
            if is_intersecting(e1, e2):
                e1.is_intersecting = True
                e2.is_intersecting = True
                

def delete_intersecting_edges(edges: List[Edge]):
    update_edges_coordinate_properties(edges)
    new_edges:List[Edge] = []
    num_edges = len(edges)
    for i in range(num_edges):
        e1 = edges[i]
        for j in range(i + 1, num_edges):
            e2 = edges[j]
            if is_intersecting(e1, e2):
                new_edges.append(e1)
                break
    return new_edges


def render(self, screen: pygame.Surface, radius: float, colour: tuple):
    screen_position = world_to_screen(self.position, WORLD_SIZE, MAGNIFICATION)
    pygame.draw.circle(screen, colour, (int(screen_position.x), int(screen_position.y)), radius)
    
    
    
def render_solid(solid: Solid,  dxf_edges: list[Edge], screen: pygame.Surface):
    for edge in solid.edges:
        edge.render(screen=screen, colour=EDGE_COLOUR_SELECTED if edge.is_intersecting else EDGE_COLOUR, width=EDGE_WIDTH)
        edge.is_intersecting = False
    for node in solid.nodes:
        node.render(screen=screen, radius=NODE_RADIUS, colour=NODE_COLOUR)
    if RENDER_DXF:
        for edge in dxf_edges:
            edge.render(screen=screen, colour=EDGE_COLOUR_SELECTED, width=EDGE_WIDTH)

def render_solid_low_level(solid:Solid, dxf_edges: list[Edge], screen: pygame.Surface):
    for i in range(solid.num_edges):
        node_1_id = solid.node_edges_connection[2*i]
        node_2_id = solid.node_edges_connection[2*i+1]
        print(f"Rendering edge {i} between nodes {node_1_id} and {node_2_id}")
        pos1 = solid.node_positions[2*node_1_id:2*node_1_id+2]
        pos2 = solid.node_positions[2*node_2_id:2*node_2_id+2]
        node1_pos = Vector(pos1[0], pos1[1])
        node2_pos = Vector(pos2[0], pos2[1])
        node1_screen_pos = world_to_screen(node1_pos, WORLD_SIZE, MAGNIFICATION)
        node2_screen_pos = world_to_screen(node2_pos, WORLD_SIZE, MAGNIFICATION)
        pygame.draw.line(screen, EDGE_COLOUR_SELECTED if edge.is_intersecting else EDGE_COLOUR, (int(node1_screen_pos.x), int(node1_screen_pos.y)), (int(node2_screen_pos.x), int(node2_screen_pos.y)), width=EDGE_WIDTH)
        edge.is_intersecting = False
    for i in range(solid.num_nodes):
        pos = solid.node_positions[2*i:2*i+2]
        screen_pos = world_to_screen(Vector(pos[0], pos[1]), WORLD_SIZE, MAGNIFICATION)
        pygame.draw.circle(screen, NODE_COLOUR, (int(screen_pos.x), int(screen_pos.y)), NODE_RADIUS)
    if RENDER_DXF:
        for edge in dxf_edges:
            edge.render(screen=screen, colour=EDGE_COLOUR_SELECTED, width=EDGE_WIDTH)

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

        
        
def load_dxf_file(file_path: str):
    doc = dxf.readfile(file_path)
    msp = doc.modelspace()
    return msp, doc

# hinge_points = []
def get_bounds_and_lines_and_hinge_points(filename:str, min_distance:float, max_distance:float, center:bool, size:float):
    msp, doc = load_dxf_file(filename)
    lines_vec:list[Tuple[Vector, Vector]] = []
    lines = []
    hinge_points:list[Vector] = []

    for e in msp.query("INSERT"):
        block = e.block()
        m = e.matrix44()
        for line in block.query("LINE"):
            start = line.dxf.start
            end = line.dxf.end
            start = m.transform(start)
            end = m.transform(end)
            lines.append((start, end))
            lines_vec.append((Vector(start.x, start.y), Vector(end.x, end.y)))
    for e in msp.query("LINE"):
        start = e.dxf.start
        end = e.dxf.end
        lines.append((start, end))
        lines_vec.append((Vector(start.x, start.y), Vector(end.x, end.y)))
    
    for i in range(len(lines_vec)):
        for j in range(i, len(lines_vec)):
            A, B, dist = nearest_points_on_line_segments(lines_vec[i][0], lines_vec[i][1], lines_vec[j][0], lines_vec[j][1])
            if min_distance < dist < max_distance:
                hinge_point = (A + B) / 2
                hinge_points.append(hinge_point)
    minx = float("inf")
    maxx = -float("inf")
    miny = float("inf")
    maxy = -float("inf")

    for i in lines_vec:
        if i[0].x < minx: minx = i[0].x
        if i[0].x > maxx: maxx = i[0].x
        if i[0].y < miny: miny = i[0].y
        if i[0].y > maxy: maxy = i[0].y

    bounds = (minx, maxx, miny, maxy)
    if center:
        cx = (minx + maxx) / 2
        cy = (miny + maxy) / 2
        bounds = (-cx, cx, -cy, cy)
        for hp in hinge_points:
            hp.x -= cx
            hp.y -= cy
        for i in lines_vec:
            i[0].x -= cx
            i[0].y -= cy
            i[1].x -= cx
            i[1].y -= cy
    if size:
        scaling = size / (maxx - minx)
        for hp in hinge_points:
            hp.x *= scaling
            hp.y *= scaling
        for i in lines_vec:
            i[0].x *= scaling
            i[0].y *= scaling
            i[1].x *= scaling
            i[1].y *= scaling
        bounds = (bounds[0] * scaling, bounds[1] * scaling, bounds[2] * scaling, bounds[3] * scaling)
    return bounds, lines_vec, hinge_points


def get_constraint_indices_and_constraint_functions(nodes:list[Node], bounds: Tuple[float, float, float, float]) -> list[Callable[[float], float]]:
    constraint_functions = []
    constraint_indices = []
    for node in nodes:
        left = (node.position.x - (bounds[0] + 2*EPSILON_BIG))<0
        right = (node.position.x - (bounds[1] - 2*EPSILON_BIG))>0
        up=(node.position.y - (bounds[2] + 2*EPSILON_BIG))<0
        down=(node.position.y - (bounds[3] - 2*EPSILON_BIG))>0
        NONBOUNDED_MIDDLE_SIZE = 0.2
        middle = np.abs(node.position.y - (bounds[2] + bounds[3]) / 2) < NONBOUNDED_MIDDLE_SIZE*(bounds[3]-bounds[2])
        if middle or up or down:
            continue
        if left:
            constraint_functions.append(lambda t: node.position.x - t*PARAMETER_MULTIPLIER)
            constraint_indices.append(2*node.id)
        elif right:
            constraint_functions.append(lambda t: node.position.x + t*PARAMETER_MULTIPLIER)
            constraint_indices.append(2*node.id + 1)
        # else:
        #     raise RuntimeError(f"Unreachable")
    return constraint_indices, constraint_functions