import pygame
from Objects import Vector, Node
SCREEN_X = 800
SCREEN_Y = 600

NODE_RADIUS = 5
NODE_COLOUR = (255, 255, 255)

def world_to_screen(position: Vector) -> Vector:
    return Vector(position.x, SCREEN_Y - position.y)
if __name__ == "__main__":
    pygame.init()
    screen = pygame.display.set_mode((SCREEN_X, SCREEN_Y))
    
    nodes = [Node(Vector(x, y)) for x in range(0, SCREEN_X, 100) for y in range(0, SCREEN_Y, 100)]
    
    
    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        screen.fill((0, 0, 0))
        for node in nodes:
            node.render(screen=screen, radius=NODE_RADIUS, colour=NODE_COLOUR, world_to_screen_function=world_to_screen)
        pygame.display.flip()

    pygame.quit()

