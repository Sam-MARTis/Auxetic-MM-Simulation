import pygame
from Objects import Vector, Node
from constants import *
from Aux import initialize_nodes



if __name__ == "__main__":
    pygame.init()
    screen = pygame.display.set_mode((screen_size, screen_size))

    nodes = initialize_nodes(100)
    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        screen.fill((0, 0, 0))
        for node in nodes:
            node.render(screen=screen, radius=NODE_RADIUS, colour=NODE_COLOUR)
        pygame.display.flip()

    pygame.quit()

