import pygame
from Objects import Vector, Node
from constants import *
from Aux import initialize_edges_distance_based, initialize_nodes, step_nodes, apply_edge_forces, apply_gravity



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
        step_nodes(nodes, delta_time)
        apply_edge_forces(edges, dt=delta_time)
        apply_gravity(nodes, dt=delta_time)

        screen.fill((0, 0, 0))
        for node in nodes:
            node.render(screen=screen, radius=NODE_RADIUS, colour=NODE_COLOUR)
        for edge in edges:
            edge.render(screen=screen)
        pygame.display.flip()

    pygame.quit()

