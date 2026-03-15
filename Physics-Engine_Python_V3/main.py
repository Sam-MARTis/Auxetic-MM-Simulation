import pygame
import time
from Objects import *
from constants import *
from Aux import *
import numba as nb



def do_animation(solid: Solid, dxf_edges:list[np.ndarray]):
    pygame.init()
    screen = pygame.display.set_mode((screen_size, screen_size))
    pygame.display.set_caption("Spring Network")
    plot_screen = pygame.Surface((PLOT_W, PLOT_H))

    t_hist: list[float] = []
    e_hist: list[float] = []
    running = True
    clock = pygame.time.Clock()
    i = 0
    stop_auto_sim = False
    @profile_function
    def anim_loop():
        nonlocal running, i, clock, stop_auto_sim
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            if event.type == pygame.KEYDOWN:
    
                left_right = False
                if event.key == pygame.K_LEFT:
                    solid.t -= DT
                    left_right = True
                    t_hist.pop()
                    e_hist.pop()
                    
                    
                if event.key == pygame.K_RIGHT:
                    solid.t += DT
                    left_right = True
                    t_hist.append(0)
                    e_hist.append(0)
                    
                if left_right:
                    # DT = DT
                    energy = solid.step(substeps=SUBSTEPS, report_energy=True)
                    t_hist[-1] = solid.t
                    e_hist[-1] = energy
                    print(f"t={solid.t:.4f}  E={energy:.6e}")
                
                

        if(solid.t > MAX_T):
            stop_auto_sim = True

        if not stop_auto_sim:

            energy = solid.step(substeps=SUBSTEPS, report_energy=True)
            solid.t += DT
        
            t_hist.append(solid.t)
            e_hist.append(energy)
            print(f"t={solid.t:.4f}  E={energy:.6e}")

        screen.fill((0, 0, 0))
        render_solid_low_level(solid=solid, dxf_edges=dxf_edges, screen=screen)
        

        
        draw_energy_plot(plot_screen, t_hist, e_hist)
        screen.blit(plot_screen, (screen_size - PLOT_W - 5, 5))
        

        
        pygame.display.flip()
    while running:
        anim_loop()
        

        
    pygame.quit()
    
def perform_sim(solid:Solid, t_min:float, t_max:float, tracking_points: np.ndarray, render: bool, dxf_edges:list[np.ndarray]):


    solid.t = t_min
    t_hist = []
    e_hist = []
    running = True
    node_positions = np.reshape(solid.node_positions, (-1, 2))
    # tracking_points_positions = [[(a[0], a[1])] for a in node_positions[tracking_points]]
    if render:
        #Set up he rendering stuff. Pygame init and all that.
        pygame.init()
        screen = pygame.display.set_mode((screen_size, screen_size))
        pygame.display.set_caption("Spring Network")
        plot_screen = pygame.Surface((PLOT_W, PLOT_H))
    while(solid.t<t_max):
        energy = solid.step(substeps=SUBSTEPS, report_energy=True)
        t_hist.append(solid.t)
        e_hist.append(energy)
        # tracked_positions = node_positions[tracking_points]
        # tracking_points_positions.append(tracked_positions)
        # for i in range(len(tracking_points)):
        #     tracking_points_positions[i].append((tracked_positions[i][0], tracked_positions[i][1]))
        solid.log_charged_nodes_positions()
        if render:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
            
                    
            # mapped_pos_points = world_to_screen_02(tracked_positions, world_size=WORLD_SIZE, magnification=MAGNIFICATION)
            # tracking_points_rendering_pos.append((mapped_pos_points[0], mapped_pos_points[1]))
            screen.fill((0, 0, 0))
            render_solid_low_level(solid=solid, dxf_edges=dxf_edges, screen=screen)
            draw_energy_plot(plot_screen=plot_screen, t_hist=t_hist, e_hist= e_hist)
            # for i in range(len(tracking_points)):
            #     render_tracking_points_pos = world_to_screen_03(tracking_points_positions[i], world_size=WORLD_SIZE, magnification=MAGNIFICATION)
            #     pygame.draw.lines(surface=screen, color=TRACKING_COLOUR, points=render_tracking_points_pos, closed=False, width=TRACKING_EDGE_WIDTH)
            solid.render_charged_points_trajectory(screen)
            screen.blit(plot_screen, (screen_size - PLOT_W - 5, 5))
            pygame.display.flip()
            if not running:
                pygame.quit()
                break
        solid.t+= DT
        


if __name__ == "__main__":
    nb.set_num_threads(17)
    min_distance = 1e-6
    max_distance = 5e-2
    center = True
    size = WORLD_SIZE/2
    n=2
    file_pattern_2cross2 = f"/home/odin/Academic/Research/Metamaterials/Auxetic Mechanical MM/Bistable/Computational Exploration/Modules/pattern_{n}x{n}.dxf"
    bounds, lines_vec, hinge_points = get_bounds_and_lines_and_hinge_points(filename=file_pattern_2cross2, min_distance=min_distance, max_distance=max_distance, center=center, size=size)
    dxf_edges = []
    for line in lines_vec:
        a = line[0]
        b = line[1]
        dxf_edges.append([np.array([a.x, a.y]), np.array([b.x, b.y])])
        
    
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

    solid.initialize_edges_distance_based(distance_threshold=LINK_DISTANCE, edges_external=dxf_edges)

    
    constraint_indices, constraint_functions = get_constraint_indices_and_constraint_functions(nodes=solid.nodes, bounds=bounds)
    solid.set_constraint(coordinate_indices_to_constrain=constraint_indices, curve_functions=constraint_functions)

    
    
    # do_animation(solid, dxf_edges=dxf_edges)
    tracking_points = np.array([np.random.randint(0, solid.num_nodes) for _ in range(500)], dtype=int)
    solid.assign_charged_indices(tracking_points)
    solid.log_charged_nodes_positions()
    perform_sim(solid=solid, t_min=0, t_max=MAX_T, tracking_points=tracking_points, render=True, dxf_edges=dxf_edges)


