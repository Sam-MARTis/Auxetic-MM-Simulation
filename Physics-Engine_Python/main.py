import pygame
from Objects import *
from constants import *
from Aux import *



if __name__ == "__main__":
    pygame.init()
    screen = pygame.display.set_mode((screen_size, screen_size))
    
    nodes = initialize_nodes(NODE_COUNT)
    edges = initialize_edges_distance_based(NODE_COUNT, nodes)
    running = True
    clock = pygame.time.Clock()
    delta_time = 0
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        delta_time = clock.tick(60) / 1000  # seconds passed since last frame
        substeps_dt = delta_time / SUBSTEPS
        for _ in range(SUBSTEPS):
            step_nodes(nodes, substeps_dt)
            apply_edge_forces(edges, dt=substeps_dt)
            apply_gravity(nodes, dt=substeps_dt)
        
        screen.fill((0, 0, 0))

        all_edges_intersection_check(edges)
        render_everything(nodes, edges, screen)

        pygame.display.flip()

    pygame.quit()

