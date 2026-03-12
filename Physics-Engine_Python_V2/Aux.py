from Objects import Vector, Node, Helper, Solid, world_to_screen
from typing import Callable, List, Tuple
from constants import *
import numpy as np
import pygame
import ezdxf as dxf
import time 
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

@profile_function
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



def draw_energy_plot(plot_screen: pygame.Surface, t_hist: list, e_hist: list):
    """Redraw the energy-vs-t plot on plot_screen."""
    plot_screen.fill((15, 15, 30))
    W, H = plot_screen.get_size()
    iw = W - 2*PLOT_PAD   
    ih = H - 2*PLOT_PAD   
    pygame.draw.line(plot_screen, (180,180,180), (PLOT_PAD, PLOT_PAD), (PLOT_PAD, PLOT_PAD+ih), 1)
    pygame.draw.line(plot_screen, (180,180,180), (PLOT_PAD, PLOT_PAD+ih), (PLOT_PAD+iw, PLOT_PAD+ih), 1)

    if len(t_hist) < 2:
        return

    t_min, t_max = t_hist[0],  max(t_hist[-1], t_hist[0]+1e-9)
    e_min, e_max = min(e_hist), max(e_hist)
    if e_max == e_min:
        e_max = e_min + 1.0   

    def to_px(t, e):
        x = PLOT_PAD + int((t - t_min) / (t_max - t_min) * iw)
        y = PLOT_PAD + ih - int((e - e_min) / (e_max - e_min) * ih)
        return x, y

    pts = [to_px(t, e) for t, e in zip(t_hist, e_hist)]
    if len(pts) >= 2:
        pygame.draw.lines(plot_screen, (0, 220, 180), False, pts, 2)

    font = pygame.font.SysFont(None, 18)
    plot_screen.blit(font.render(f"t={t_hist[-1]:.3f}", True, (200,200,200)), (PLOT_PAD+iw-60, PLOT_PAD+ih+4))
    plot_screen.blit(font.render(f"E={e_hist[-1]:.3e}", True, (200,200,200)), (2, PLOT_PAD))
    plot_screen.blit(font.render("Energy vs t", True, (220,220,220)), (W//2-30, 4))



def render(self, screen: pygame.Surface, radius: float, colour: tuple):
    screen_position = world_to_screen(self.position, WORLD_SIZE, MAGNIFICATION)
    pygame.draw.circle(screen, colour, (int(screen_position.x), int(screen_position.y)), radius)
    
    
    

def map_strain_to_pygame_colour(strain: float) -> tuple:
    strain_scaled = 50000*strain
    r = max(min(int(125+strain_scaled), 255), 0)
    b = max(min(int(125-strain_scaled), 255), 0)
    g = 0
    return (r, g, b)

def render_solid_low_level(solid:Solid, dxf_edges: list[np.ndarray], screen: pygame.Surface):
    for i in range(solid.num_edges):
        
        node_1_id = solid.node_edges_connection[2*i]
        node_2_id = solid.node_edges_connection[2*i+1]
        
        pos1 = solid.node_positions[2*node_1_id:2*node_1_id+2]
        pos2 = solid.node_positions[2*node_2_id:2*node_2_id+2]
        node1_pos = Vector(pos1[0], pos1[1])
        node2_pos = Vector(pos2[0], pos2[1])
        node1_screen_pos = world_to_screen(node1_pos, WORLD_SIZE, MAGNIFICATION)
        node2_screen_pos = world_to_screen(node2_pos, WORLD_SIZE, MAGNIFICATION)

        l0 = solid.edges_extension_natural[i]
        colour = map_strain_to_pygame_colour((solid.edges_extension[i]-l0)/l0)

        pygame.draw.line(screen, colour, (int(node1_screen_pos.x), int(node1_screen_pos.y)), (int(node2_screen_pos.x), int(node2_screen_pos.y)), width=EDGE_WIDTH)
    for i in range(solid.num_nodes):
        
        pos = solid.node_positions[2*i:2*i+2]
        constrained = (2*i in solid.fixed_indices) or (2*i+1 in solid.fixed_indices) 
        screen_pos = world_to_screen(Vector(pos[0], pos[1]), WORLD_SIZE, MAGNIFICATION)
        pygame.draw.circle(screen, NODE_COLOUR if not constrained else NODE_COLOUR_SELECTED, (int(screen_pos.x), int(screen_pos.y)), NODE_RADIUS)
    if RENDER_DXF:
        for i in range(len(dxf_edges)):
            a = dxf_edges[i][0]
            b = dxf_edges[i][1]
            a_screen = world_to_screen(Vector(a[0], a[1]), WORLD_SIZE, MAGNIFICATION)
            b_screen = world_to_screen(Vector(b[0], b[1]), WORLD_SIZE, MAGNIFICATION)
            pygame.draw.line(screen, EDGE_COLOUR_SELECTED, (int(a_screen.x), int(a_screen.y)), (int(b_screen.x), int(b_screen.y)), width=1)

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

@profile_function
def get_constraint_indices_and_constraint_functions(nodes:list[Node], bounds: Tuple[float, float, float, float]) -> list[Callable[[float], float]]:
    constraint_functions = []
    constraint_indices = []
    for node in nodes:
        left  = node.position.x < bounds[0] + 2*EPSILON_BIG
        right = node.position.x > bounds[1] - 2*EPSILON_BIG
        down = node.position.y < bounds[2] + 2*EPSILON_BIG
        NONBOUNDED_MIDDLE_SIZE = 0.1

        middle = np.abs(node.position.y - (bounds[2] + bounds[3]) / 2) < NONBOUNDED_MIDDLE_SIZE*(bounds[3]-bounds[2])
        if middle:
            continue
        if left:
            x0 = node.position.x
            y0 = node.position.y
            constraint_functions.append(lambda t, _x=x0: _x - t*PARAMETER_MULTIPLIER)
            constraint_indices.append(2*node.id)
            # constraint_functions.append(lambda t, _y=y0: _y)
            # constraint_indices.append(2*node.id + 1)
        elif right:
            x0 = node.position.x
            y0 = node.position.y
            constraint_functions.append(lambda t, _x=x0: _x + t*PARAMETER_MULTIPLIER)
            constraint_indices.append(2*node.id)
            # constraint_functions.append(lambda t, _y=y0: _y)
            # constraint_indices.append(2*node.id + 1)
        elif down:
            y0 = node.position.y
            constraint_functions.append(lambda t, _y=y0: _y)
            constraint_indices.append(2*node.id + 1)
    return constraint_indices, constraint_functions