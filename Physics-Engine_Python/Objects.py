from typing import List, Callable
import pygame


class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def dot(self, other: 'Vector') -> float:
        return self.x * other.x + self.y * other.y
    def magnitude(self) -> float:
        return (self.x ** 2 + self.y ** 2) ** 0.5

class Node:
    def __init__(self, position: Vector):
        self.connections = []
        self.position = position
    def render(self, screen: pygame.Surface, radius: float, colour: tuple, world_to_screen_function: Callable[[Vector], Vector]):
        screen_position = world_to_screen_function(self.position)
        pygame.draw.circle(screen, colour, (int(screen_position.x), int(screen_position.y)), radius)

