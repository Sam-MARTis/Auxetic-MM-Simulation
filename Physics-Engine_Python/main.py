import pygame
from Objects import *
from constants import *
from Aux import *



def do_animation(nodes: list[Node], edges: list[Edge]):
    
    pygame.init()
    gravity = GRAVITY
    hold = True
    screen = pygame.display.set_mode((screen_size, screen_size))
    
    node_selected_by_mouse = None
    running = True
    clock = pygame.time.Clock()
    delta_time = 0
    i = 0
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
            # Mouseclick once
            if event.type == pygame.MOUSEBUTTONDOWN:
                mouse_pos = pygame.mouse.get_pos()
                for node in nodes:
                    node_screen_pos = world_to_screen(node.position, WORLD_SIZE, MAGNIFICATION)
                    dx = mouse_pos[0] - node_screen_pos.x
                    dy = mouse_pos[1] - node_screen_pos.y
                    if dx**2 + dy**2 < 2*NODE_RADIUS**2:
                        node_selected_by_mouse = node
                        break
                print("MOUSEDOWN")
                if node_selected_by_mouse:
                    print(f"Node selected: {node_selected_by_mouse}")
            if event.type == pygame.MOUSEBUTTONUP:
                if node_selected_by_mouse:
                    print(f"Node released: {node_selected_by_mouse}")
                node_selected_by_mouse = None
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_1:
                    edges[i].delete_edge()
                    i+=1
                    print("Deleted edge")

        delta_time = clock.tick(60) / 1000  # seconds passed since last frame
        substeps_dt = delta_time / SUBSTEPS
        for _ in range(SUBSTEPS):
            
            step_nodes(nodes, substeps_dt)
            if hold:
                for node in nodes[-5:]:
                    node.position.y =0.5
            if node_selected_by_mouse:
                mouse_pos = pygame.mouse.get_pos()
                world_x = mouse_pos[0] / MAGNIFICATION - WORLD_SIZE / 2
                world_y = -(mouse_pos[1] / MAGNIFICATION - WORLD_SIZE / 2)
                node_selected_by_mouse.position.x = world_x
                node_selected_by_mouse.position.y = world_y
            apply_edge_forces(edges, dt=substeps_dt)
            apply_gravity(nodes, dt=substeps_dt, gravity=gravity)
        
        screen.fill((0, 0, 0))

        render_everything(nodes, edges, screen)
        

        pygame.display.flip()

    pygame.quit()

if __name__ == "__main__":

    nodes = initialize_nodes(NODE_COUNT)
    nodes[-5:] = [Node(mass=NODE_MASS, position=Vector(i/10, 0.5), velocity=Vector(0, 0)) for i in range(5)]
    edges = initialize_edges_distance_based(NODE_COUNT, nodes)
    
    # edges = all_edges_intersection_check(edges)
    do_animation(nodes, edges)
    

