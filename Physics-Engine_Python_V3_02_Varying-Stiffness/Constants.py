import numpy as np

PROFILE_PERFORMANCE = False

EPSILON = 1e-6
EPSILON_BIG = 1e-4
TOLERANCE_SOLVER = 1e-1
MAX_ITER_SOLVER = 280
RELAXATION = 1.0
RELAXATION2 = 0.1

DT = 0.5e-0
PARAMETER_MULTIPLIER = 3e-2
MAX_T = 25
ALPHA_0 = 1e-2
DIMS = 2

CHARGE_OPTIMIZATION_FILE = "charge_optimize.py"
DISTANCE_DATA_STORAGE_FILE = "charge_distance_data.csv"
OPTIMIZED_CHARGE_STORAGE_FILE = "optimized_charges.csv" 
OPTIMIZATION_VERBOSE = True

OBJECT_SAVE_FILE = "solid_state.pkl"


PLOT_W, PLOT_H = 500, 300
PLOT_PAD = 40   



WORLD_SIZE = 10
MAGNIFICATION = 100
NODE_RADIUS = 1
NODE_COLOUR = (255, 255, 255)
NODE_COLOUR_SELECTED = (255, 200, 0)
TRACKING_COLOUR = (100, 255, 0)
NODE_MASS = 1
NODE_COUNT = 2000
DENSITY_FACTOR = 3.5
STARTING_VEL = 0.5

DXF_EDGE_WIDTH = 1
TRACKING_EDGE_WIDTH = 2

SCALING = 1
VISCOCITY = 10*SCALING

SHAPE_SIZE = WORLD_SIZE/2

GRAVITY_ENABLED = False
GRAVITY = 1.81

CHARGE_CONSTANT = 0.1
SUBSTEPS = 1

BOUNDARY_X_P = WORLD_SIZE/2.1
BOUNDARY_X_N = -WORLD_SIZE/2.1
BOUNDARY_Y_P = WORLD_SIZE/2.1
BOUNDARY_Y_N = -WORLD_SIZE/2.1

adjacent_bond_size = SHAPE_SIZE/(np.sqrt(NODE_COUNT)-1)
# For distance based attachment of links
LINK_DISTANCE =  adjacent_bond_size*np.sqrt(2.2)

EDGE_COLOUR = (0, 200, 200)
EDGE_COLOUR_SELECTED = (200, 100, 0)
EDGE_WIDTH = 2
LINK_STRENGTH = 5000*SCALING

screen_size = WORLD_SIZE * (MAGNIFICATION)
RENDER_DXF = True

NUM_NON_CHARGED_PER_CHARGED = 10





import time

def profile_function(func):
    if PROFILE_PERFORMANCE:
        def wrapper(*args, **kwargs):
            start_time = time.perf_counter()
            result = func(*args, **kwargs)
            end_time = time.perf_counter()
            print(f"{func.__name__} took {(end_time - start_time) * 1e3:.2f} ms")
            return result
        return wrapper
    return func
