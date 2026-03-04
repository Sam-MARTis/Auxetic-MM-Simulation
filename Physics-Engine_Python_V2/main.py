import pygame
from Objects import *
from constants import *
from Aux import *
import numba as nb



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
                solid.step(dt=delta_time, substeps=SUBSTEPS)
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

        delta_time = DT  # seconds passed since last frame

        # solid.step(dt=delta_time, substeps=SUBSTEPS)

        if node_selected_by_mouse:
            mouse_pos = pygame.mouse.get_pos()
            world_x = mouse_pos[0] / MAGNIFICATION - WORLD_SIZE / 2
            world_y = -(mouse_pos[1] / MAGNIFICATION - WORLD_SIZE / 2)
            solid.node_positions[2*i] = world_x
            solid.node_positions[2*i+1] = world_y

        screen.fill((0, 0, 0))

        render_solid_low_level(solid=solid, dxf_edges=dxf_edges, screen=screen)

        pygame.display.flip()

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
    

