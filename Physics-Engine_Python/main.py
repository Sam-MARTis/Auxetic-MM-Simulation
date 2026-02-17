import pygame
from Objects import *
from constants import *
from Aux import *



def do_animation(solid: Solid, dxf_edges: list[Edge]):
    pygame.init()
    gravity = GRAVITY
    hold = True
    screen = pygame.display.set_mode((screen_size, screen_size))
    
    node_selected_by_mouse = None
    del_edge_key_clicked = False
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
                for node in solid.nodes:
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

            if del_edge_key_clicked:
                solid.delete_edge(solid.edges[0])
                # i+=1
                print("Deleted edge")
                print(f"Remaining edges: {len(solid.edges)}")

        delta_time = clock.tick(60) / 1000  # seconds passed since last frame

        solid.step(dt=delta_time, substeps=SUBSTEPS, gravity=gravity)

        if node_selected_by_mouse:
            mouse_pos = pygame.mouse.get_pos()
            world_x = mouse_pos[0] / MAGNIFICATION - WORLD_SIZE / 2
            world_y = -(mouse_pos[1] / MAGNIFICATION - WORLD_SIZE / 2)
            node_selected_by_mouse.position.x = world_x
            node_selected_by_mouse.position.y = world_y
        
        screen.fill((0, 0, 0))

        render_solid(solid=solid, dxf_edges=dxf_edges, screen=screen)

        pygame.display.flip()

    pygame.quit()

if __name__ == "__main__":
    min_distance = 1e-6
    max_distance = 5e-2
    center = True
    size = WORLD_SIZE/2
    bounds, lines_vec, hinge_points = get_bounds_and_lines_and_hinge_points(filename="/home/odin/Academic/Research/Metamaterials/Auxetic Mechanical MM/Bistable/Computational Exploration/Modules/mm_pattern_01.dxf", min_distance=min_distance, max_distance=max_distance, center=center, size=size)
    dxf_edges:list[Edge] = []
    for line in lines_vec:
        a = line[0]
        b = line[1]
        print(f"Line from {a} to {b}")
        node_a = Node(mass=NODE_MASS, position=Vector(a.x, a.y), velocity=Vector(0, 0))
        node_b = Node(mass=NODE_MASS, position=Vector(b.x, b.y), velocity=Vector(0, 0))
        edge = Edge(k=LINK_STRENGTH, l0=(b - a).magnitude(), node_a=node_a, node_b=node_b, attach=True)
        dxf_edges.append(edge)
        print(node_b.position)
    hinge_point_nodes = [Node(mass=NODE_MASS, position=hp, velocity=Vector(0, 0)) for hp in hinge_points]
    nodes = initialize_nodes_uniform_with_noise(bounds=bounds, node_count_x=int(np.sqrt(NODE_COUNT)), node_count_y=int(np.sqrt(NODE_COUNT)))
    nodes.extend(hinge_point_nodes)
    solid = Solid(nodes=nodes)
    
    solid.initialize_edges_distance_based(distance_threshold=LINK_DISTANCE)

    solid.remove_intersecting_edges(dxf_edges)
    solid.remove_intersecting_edges(dxf_edges)
    
    do_animation(solid, dxf_edges=dxf_edges)
    

