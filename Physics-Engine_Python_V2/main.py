import pygame
import time
from Objects import *
from constants import *
from Aux import *
import numba as nb



def draw_energy_plot(plot_screen: pygame.Surface, t_hist: list, e_hist: list):
    """Redraw the energy-vs-t plot on plot_screen."""
    plot_screen.fill((15, 15, 30))
    W, H = plot_screen.get_size()
    iw = W - 2*PLOT_PAD   # inner width
    ih = H - 2*PLOT_PAD   # inner height

    # Axes
    pygame.draw.line(plot_screen, (180,180,180), (PLOT_PAD, PLOT_PAD), (PLOT_PAD, PLOT_PAD+ih), 1)
    pygame.draw.line(plot_screen, (180,180,180), (PLOT_PAD, PLOT_PAD+ih), (PLOT_PAD+iw, PLOT_PAD+ih), 1)

    if len(t_hist) < 2:
        return

    t_min, t_max = t_hist[0],  max(t_hist[-1], t_hist[0]+1e-9)
    e_min, e_max = min(e_hist), max(e_hist)
    if e_max == e_min:
        e_max = e_min + 1.0   # avoid degenerate range

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


def do_animation(solid: Solid, dxf_edges: list[Edge]):
    pygame.init()
    gravity = GRAVITY
    hold = True
    screen = pygame.display.set_mode((screen_size, screen_size))
    pygame.display.set_caption("Spring Network")
    plot_screen = pygame.Surface((PLOT_W, PLOT_H))

    t_hist: list[float] = []
    e_hist: list[float] = []
    
    node_selected_by_mouse = None
    del_edge_key_clicked = False
    running = True
    clock = pygame.time.Clock()
    delta_time = 0
    i = 0
    time.sleep(3)
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            if event.type == pygame.KEYDOWN:
                # Arrow keys
                if event.key == pygame.K_UP:
                    gravity += 5
                if event.key == pygame.K_DOWN:
                    gravity -= 5
                if event.key == pygame.K_r:
                    hold = False
                        
                if event.key == pygame.K_1:
                                
                    solid.delete_edge(solid.edges[i])
                    i+=1
                    print("Deleted edge")
                    print(f"Remaining edges: {len(solid.edges)}")

                if event.key == pygame.K_1:
                    del_edge_key_clicked = True
            if event.type == pygame.KEYUP:
                if event.key == pygame.K_1:
                    del_edge_key_clicked = False
            # Mouseclick once
            if event.type == pygame.MOUSEBUTTONDOWN:
                
                mouse_pos = pygame.mouse.get_pos()
                for i in range(len(solid.nodes)):
                    node_pos = Vector(solid.node_positions[2*i], solid.node_positions[2*i+1])
                    node_screen_pos = world_to_screen(node_pos, WORLD_SIZE, MAGNIFICATION)
                    dx = mouse_pos[0] - node_screen_pos.x
                    dy = mouse_pos[1] - node_screen_pos.y
                    if dx**2 + dy**2 < 2*NODE_RADIUS**2:
                        node_selected_by_mouse = i
                        break
                print("MOUSEDOWN")
                if node_selected_by_mouse:
                    print(f"Node selected: {node_selected_by_mouse}")
            if event.type == pygame.MOUSEBUTTONUP:
                if node_selected_by_mouse:
                    print(f"Node released: {node_selected_by_mouse}")
                node_selected_by_mouse = None

            if del_edge_key_clicked:
                solid.delete_edge(solid.edges[0])
                # i+=1
                print("Deleted edge")
                print(f"Remaining edges: {len(solid.edges)}")
        if node_selected_by_mouse:
            mouse_pos = pygame.mouse.get_pos()
            world_x = mouse_pos[0] / MAGNIFICATION - WORLD_SIZE / 2
            world_y = -(mouse_pos[1] / MAGNIFICATION - WORLD_SIZE / 2)
            solid.node_positions[2*i] = world_x
            solid.node_positions[2*i+1] = world_y

        if(solid.t < MAX_T):
            delta_time = DT
            _s = time.perf_counter() if PROFILE_PERFORMANCE else 0
            energy = solid.step(dt=delta_time, substeps=SUBSTEPS, report_energy=True)
            _step_ms = (time.perf_counter() - _s) * 1e3 if PROFILE_PERFORMANCE else 0
            t_hist.append(solid.t)
            e_hist.append(energy)
            print(f"t={solid.t:.4f}  E={energy:.6e}")

        _s = time.perf_counter() if PROFILE_PERFORMANCE else 0
        screen.fill((0, 0, 0))
        render_solid_low_level(solid=solid, dxf_edges=dxf_edges, screen=screen)
        _render_ms = (time.perf_counter() - _s) * 1e3 if PROFILE_PERFORMANCE else 0

        _s = time.perf_counter() if PROFILE_PERFORMANCE else 0
        draw_energy_plot(plot_screen, t_hist, e_hist)
        screen.blit(plot_screen, (screen_size - PLOT_W - 5, 5))
        _plot_ms = (time.perf_counter() - _s) * 1e3 if PROFILE_PERFORMANCE else 0

        _s = time.perf_counter() if PROFILE_PERFORMANCE else 0
        pygame.display.flip()
        _flip_ms = (time.perf_counter() - _s) * 1e3 if PROFILE_PERFORMANCE else 0

        if PROFILE_PERFORMANCE:
            print(f"  [render]  solid={_render_ms:.2f}ms  plot={_plot_ms:.2f}ms  flip={_flip_ms:.2f}ms  | step={_step_ms:.2f}ms")

    pygame.quit()

if __name__ == "__main__":
    nb.set_num_threads(17)
    min_distance = 1e-6
    max_distance = 5e-2
    center = True
    size = WORLD_SIZE/2
    n=2
    file_pattern_2cross2 = f"/home/odin/Academic/Research/Metamaterials/Auxetic Mechanical MM/Bistable/Computational Exploration/Modules/pattern_{n}x{n}.dxf"
    bounds, lines_vec, hinge_points = get_bounds_and_lines_and_hinge_points(filename=file_pattern_2cross2, min_distance=min_distance, max_distance=max_distance, center=center, size=size)
    dxf_edges:list[Edge] = []
    for line in lines_vec:
        a = line[0]
        b = line[1]
        # print(f"Line from {a} to {b}")
        node_a = Node(mass=NODE_MASS, position=Vector(a.x, a.y), velocity=Vector(0, 0))
        node_b = Node(mass=NODE_MASS, position=Vector(b.x, b.y), velocity=Vector(0, 0))
        edge = Edge(k=LINK_STRENGTH, l0=(b - a).magnitude(), node_a=node_a, node_b=node_b, attach=True)
        dxf_edges.append(edge)
        # print(node_b.position)
    
    Node._count=0
    print(f"Node class count: {Node._count}")
    nodes, boundary_nodes = initialize_nodes_uniform_with_noise(bounds=bounds, node_count_x=int(np.sqrt(NODE_COUNT)), node_count_y=int(np.sqrt(NODE_COUNT)), get_boundary_nodes=True)
    print(f"Nodes without hinge points len: {len(nodes)} vs Class node count: {Node._count}")
    hinge_point_nodes = [Node(mass=NODE_MASS, position=hp, velocity=Vector(0, 0)) for hp in hinge_points]
    print(f"hinge points len: {len(hinge_point_nodes)} vs Class node count: {Node._count}")
    nodes.extend(hinge_point_nodes)
    print(f"Node count: {len(nodes)}")
    print(f"Class node count: {Node._count}")
    solid = Solid(nodes=nodes)
    
    solid.initialize_edges_distance_based(distance_threshold=LINK_DISTANCE)

    solid.remove_intersecting_edges(dxf_edges)
    solid.remove_intersecting_edges(dxf_edges)
    constraint_indices, constraint_functions = get_constraint_indices_and_constraint_functions(nodes=solid.nodes, bounds=bounds)
    solid.prime()
    solid.set_constraint(coordinate_indices_to_constrain=constraint_indices, curve_functions=constraint_functions)
    # solid.step(dt=DT, substeps=SUBSTEPS)
    
    do_animation(solid, dxf_edges=dxf_edges)
    

