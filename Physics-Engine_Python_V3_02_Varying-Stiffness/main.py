import pygame
import time
from Objects import *
from Constants import *
from Aux import *


np.random.seed(42)

class Button:
    def __init__(self, x: float, y: float, width: float, height: float, text: str, color: tuple = (100, 100, 100), text_color: tuple = (255, 255, 255)):
        self.rect = pygame.Rect(x, y, width, height)
        self.text = text
        self.color = color
        self.text_color = text_color
        self.hover = False
        
    def draw(self, screen: pygame.Surface, font: pygame.font.Font):
        color = tuple(min(c + 50, 255) for c in self.color) if self.hover else self.color
        pygame.draw.rect(screen, color, self.rect)
        pygame.draw.rect(screen, (200, 200, 200), self.rect, 2)
        text_surface = font.render(self.text, True, self.text_color)
        text_rect = text_surface.get_rect(center=self.rect.center)
        screen.blit(text_surface, text_rect)
        
    def is_clicked(self, pos: tuple) -> bool:
        return self.rect.collidepoint(pos)
        
    def update_hover(self, pos: tuple):
        self.hover = self.rect.collidepoint(pos) 

def do_animation(solid: Solid, dxf_edges:list[np.ndarray], relaxation: float):
    pygame.init()
    screen = pygame.display.set_mode((screen_size, screen_size))
    pygame.display.set_caption("Spring Network")
    plot_screen = pygame.Surface((PLOT_W, PLOT_H))
    
    font_large = pygame.font.SysFont(None, 32)
    font_small = pygame.font.SysFont(None, 20)

    button_width, button_height = 70, 30
    button_x = 10
    button_y = 10
    pause_btn = Button(button_x, button_y, button_width, button_height, "Pause", color=(100, 100, 150))
    resume_btn = Button(button_x + 80, button_y, button_width, button_height, "Resume", color=(100, 150, 100))
    restart_btn = Button(button_x + 160, button_y, button_width, button_height, "Restart", color=(150, 100, 100))
    rewind_btn = Button(button_x + 240, button_y, button_width, button_height, "Rewind", color=(100, 100, 100))

    t_hist: list[float] = []
    e_hist: list[float] = []
    running = True
    clock = pygame.time.Clock()
    paused = False
    current_energy = 0.0
    switch_state = False
    
    @profile_function
    def anim_loop():
        nonlocal running, paused, current_energy, clock, switch_state
        
        mouse_pos = pygame.mouse.get_pos()
        pause_btn.update_hover(mouse_pos)
        resume_btn.update_hover(mouse_pos)
        restart_btn.update_hover(mouse_pos)
        rewind_btn.update_hover(mouse_pos)
        
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            
            if event.type == pygame.MOUSEBUTTONDOWN:
                if pause_btn.is_clicked(mouse_pos):
                    paused = True
                elif resume_btn.is_clicked(mouse_pos):
                    paused = False
                elif restart_btn.is_clicked(mouse_pos):
                    solid.t = 0
                    t_hist.clear()
                    e_hist.clear()
                    current_energy = 0.0
                    paused = False
                elif rewind_btn.is_clicked(mouse_pos):
                    paused = True
                    switch_state = not switch_state
            
            if event.type == pygame.KEYDOWN:
                left_right = False
                if event.key == pygame.K_LEFT:
                    if len(t_hist) > 0:
                        solid.t -= DT
                        left_right = True
                        t_hist.pop()
                        e_hist.pop()
                        if t_hist:
                            current_energy = e_hist[-1]
                    
                if event.key == pygame.K_RIGHT:
                    solid.t += DT
                    left_right = True
                    t_hist.append(0)
                    e_hist.append(0)
                    
                if left_right:
                    energy = solid.step(substeps=SUBSTEPS, report_energy=True, relaxation=relaxation)
                    t_hist[-1] = solid.t
                    e_hist[-1] = energy
                    current_energy = energy
                    print(f"t={solid.t:.4f}  E={energy:.6e}")

        if not paused and solid.t <= MAX_T:
            energy = solid.step(substeps=SUBSTEPS, report_energy=True, relaxation=relaxation)
            solid.t += DT
        
            t_hist.append(solid.t)
            e_hist.append(energy)
            current_energy = energy
            print(f"t={solid.t:.4f}  E={energy:.6e}")

        screen.fill((0, 0, 0))
        render_solid_low_level(solid=solid, dxf_edges=dxf_edges, screen=screen)

        draw_energy_plot(plot_screen, t_hist, e_hist)
        screen.blit(plot_screen, (screen_size - PLOT_W - 5, 5))
        
        pause_btn.draw(screen, font_small)
        resume_btn.draw(screen, font_small)
        restart_btn.draw(screen, font_small)
        rewind_btn.draw(screen, font_small)
        
        status_text = f"t = {solid.t:.4f}"
        energy_text = f"E = {current_energy:.6e}"
        status_surface = font_large.render(status_text, True, (0, 255, 100))
        energy_surface = font_large.render(energy_text, True, (0, 200, 255))
        
        screen.blit(status_surface, (20, 50))
        screen.blit(energy_surface, (20, 90))
        
        if paused:
            pause_text = font_large.render("PAUSED", True, (255, 100, 100))
            screen.blit(pause_text, (screen_size // 2 - 70, 20))
        
        pygame.display.flip()
        
    while running:
        anim_loop()


    # solid.t = 0
    # t_hist = []
    # e_hist = []
    # running = True
    
    
    # if render:
    #     #Set up he rendering stuff. Pygame init and all that.
    #     pygame.init()
    #     screen = pygame.display.set_mode((screen_size, screen_size))
    #     pygame.display.set_caption("Spring Network")
    #     plot_screen = pygame.Surface((PLOT_W, PLOT_H))
    # while(solid.t<MAX_T):
    #     solid.step(substeps=SUBSTEPS, report_energy=True, relaxation=relaxation)
    #     solid.log_state()
    #     if render:
    #         for event in pygame.event.get():
    #             if event.type == pygame.QUIT:
    #                 running = False
            
                    
        
    #         screen.fill((0, 0, 0))
    #         render_solid_low_level(solid=solid, dxf_edges=dxf_edges, screen=screen)
    #         draw_energy_plot(plot_screen=plot_screen, t_hist=solid.t_log, e_hist=solid.energy_log)

    #         solid.render_charged_points_trajectory(screen)
    #         screen.blit(plot_screen, (screen_size - PLOT_W - 5, 5))
    #         pygame.display.flip()
    #         if not running:
    #             pygame.quit()
    #             break
    #     solid.t+= DT
        


if __name__ == "__main__":
    # nb.set_num_threads(17)
    min_distance = 1e-6
    max_distance = 5e-2
    center = True
    size = WORLD_SIZE/2
    n=1
    file_pattern_2cross2 = f"./dxf_files/pattern_{n}x{n}_edge-gap_0.03.dxf"
    bounds, lines_vec, hinge_points = get_bounds_and_lines_and_hinge_points(filename=file_pattern_2cross2, min_distance=min_distance, max_distance=max_distance, center=center, size=size)
    dxf_edges = []
    for line in lines_vec:
        a = line[0]
        b = line[1]
        dxf_edges.append([np.array([a.x, a.y]), np.array([b.x, b.y])])
        
    
    Node._count=0
    print(f"Node class count: {Node._count}")
    nodes, boundary_nodes, charged_indices = initialize_nodes_uniform_with_noise(bounds=bounds, node_count_x=int(np.sqrt(NODE_COUNT)), node_count_y=int(np.sqrt(NODE_COUNT)), get_boundary_nodes=True)
    print(f"Nodes without hinge points len: {len(nodes)} vs Class node count: {Node._count}")
    hinge_point_nodes = [Node(mass=NODE_MASS, position=hp, velocity=Vector(0, 0)) for hp in hinge_points]
    print(f"hinge points len: {len(hinge_point_nodes)} vs Class node count: {Node._count}")
    nodes.extend(hinge_point_nodes)
    print(f"Node count: {len(nodes)}")
    print(f"Class node count: {Node._count}")
    solid = Solid(nodes=nodes)
    print(f"Bounds are: {bounds}")

    solid.initialize_edges_distance_based(distance_threshold=LINK_DISTANCE, edges_external=dxf_edges)

    
    constraint_indices, constraint_functions = get_constraint_indices_and_constraint_functions(nodes=solid.nodes, bounds=bounds)
    solid.set_constraint(coordinate_indices_to_constrain=constraint_indices, curve_functions=constraint_functions)

    
    
    # do_animation(solid, dxf_edges=dxf_edges)
    # tracking_points = np.array([np.random.randint(0, solid.num_nodes) for _ in range(50)], dtype=int)
    # solid.assign_charged_indices(charged_indices)
    # solid.log_charged_nodes_positions()


    solid.log_state()
    save_solid_to_file(solid, OBJECT_SAVE_FILE)
    # perform_sim(solid=solid, t_min=0, t_max=MAX_T, tracking_points=charged_indices, render=False, dxf_edges=dxf_edges, relaxation=RELAXATION)
    do_animation(solid, dxf_edges=dxf_edges, relaxation=RELAXATION)


    # solid.write_data_to_csv(filename=DISTANCE_DATA_STORAGE_FILE)
    # solid = load_solid_from_file(OBJECT_SAVE_FILE)
    # solid.optimize_charges(worker_file=CHARGE_OPTIMIZATION_FILE, position_storage_file=DISTANCE_DATA_STORAGE_FILE, optimized_charge_storage_file=OPTIMIZED_CHARGE_STORAGE_FILE, verbose=OPTIMIZATION_VERBOSE)
    # perform_sim(solid=solid, t_min=0, t_max=MAX_T, tracking_points=charged_indices, render=True, dxf_edges=dxf_edges, relaxation = RELAXATION2)


