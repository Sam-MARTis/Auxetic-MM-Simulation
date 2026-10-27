#pragma once

#include <limits>
#include <string>
#include <Eigen/Sparse>
#include <Eigen/Dense>
// #define NDEBUG

// Index type used for every mesh indirection (nodes, edges, boundaries).
// Matches the FVM reference so the two code bases read the same way.
using Index = unsigned int;

// Sentinel written into a tag slot that carries no boundary.
constexpr Index NO_BOUNDARY = std::numeric_limits<Index>::max();

// #define SHOULD_PROFILE

constexpr double EPSILON = 1.0e-8;
constexpr double DEFAULT_DISTANCE_FROM_CORNER = 1.0e-4;
constexpr double DEFAULT_INTERNAL_EDGE_INSERT_OFFSET = 1.0e-2;
constexpr double DEFAULT_PARTICLE_RENDER_RADIUS = 1e-2;
constexpr double DEFAULT_SCREEN_WIDTH = 800;
constexpr double DEFAULT_SCREEN_HEIGHT = 800;
constexpr double DEFAULT_SCREEN_SCALE = DEFAULT_SCREEN_WIDTH/3.5;
constexpr unsigned int DEFAULT_NUMBER_INSERT_TRIES = 100;
constexpr unsigned int DEFAULT_NUMBER_INTERNAL_NODES = 300;
constexpr double TARGET_BOUNDARY_NODE_DENSITY = 20.0;
constexpr double TARGET_INTERNAL_EDGE_NODE_DENSITY = 10.0;
constexpr unsigned int DEFAULT_NUMBER_AREA_SAMPLES = 1000;
constexpr unsigned int DEFAULT_NUMBER_NODES_PER_TREE_LEAF = 10;
constexpr unsigned int DEFAULT_NUM_NEIGHBOURS_FOR_EDGES = 15;
constexpr double DEFAULT_DISTANCE_THRESHOLD_FOR_EDGES = 12e-2;
constexpr double DEFAULT_EDGE_RENDER_THICKNESS = 4e-3;
constexpr double DEFAULT_EDGE_STIFFNESS = 1.0;
// Bond-angle (torsional) spring. Dimensionless prefactor for the per-element
// stiffness kappa = DEFAULT_TORSIONAL_STIFFNESS * 0.5*(k1+k2)*l1*l2; see
// physicsDerivation.md section 13.
constexpr double DEFAULT_TORSIONAL_STIFFNESS = 1.0;
// An element whose sine is below this is (nearly) collinear. The cosine energy
// is stationary there, so count it as a diagnostic rather than silently skipping.
constexpr double NEAR_COLLINEAR_SINE_THRESHOLD = 1e-6;
constexpr double DEFAULT_BOUNDARY_NORMAL_RENDER_LENGTH = 4e-2;
// Incident boundary chords whose normals differ by less than this dot product are
// treated as one plane; a sharper corner gets a constraint row per edge.
constexpr double BOUNDARY_CORNER_DOT_THRESHOLD = 0.966;   // cos(15 deg)

// Gmsh meshing (see include/gmsh_mesher.hpp)
constexpr bool DEFAULT_USE_GMSH_MESHER = true;
constexpr double DEFAULT_GMSH_MESH_SIZE = 1.0 / TARGET_INTERNAL_EDGE_NODE_DENSITY;
// Cuts are removed as slots of this width. Keep it near the mesh size: a thinner
// slot leaves stiff short springs at the tip and destabilises the solve.
constexpr double DEFAULT_GMSH_SLOT_WIDTH = DEFAULT_GMSH_MESH_SIZE;
constexpr bool DEFAULT_GMSH_WRITE_MESH = false;
constexpr const char* DEFAULT_GMSH_MESH_PATH = "gmsh_mesh.msh";

static constexpr int    MAX_SOLVER_INTERNAL_ITER = 5;
static constexpr double TOL      = 1e-2;
static constexpr double EPSILON2  = 1e-12;

// Energy-decreasing backtracking line search on the Newton step; see
// physicsDerivation.md section 14. The constraint is linear, so any scaled step
// keeps A dx = -violation and never increases the violation.
constexpr int    LINE_SEARCH_MAX_ITERATIONS = 20;
constexpr double LINE_SEARCH_INITIAL_STEP = 1.0;
constexpr double LINE_SEARCH_REDUCTION = 0.5;
constexpr double LINE_SEARCH_ARMIJO = 1e-4;
constexpr double LINE_SEARCH_MIN_STEP = 1e-8;
constexpr double STEP_TOL = 1e-8;
// Used only when the Hessian is indefinite and the Newton step is not a descent
// direction: K + shift*I is retried with a growing shift.
constexpr double HESSIAN_SHIFT_INITIAL = 1e-8;
constexpr double HESSIAN_SHIFT_FACTOR = 10.0;
constexpr double HESSIAN_SHIFT_MAX = 1e8;
// Above this plane violation the KKT step is a feasibility step for a moving
// constraint target; the energy may legitimately increase, so no energy line
// search is applied. After a good KKT solve the violation is ~1e-16.
constexpr double CONSTRAINT_FEASIBILITY_TOL = 1e-8;

constexpr double DEFAULT_INTERNAL_EDGE_EXPANSIVE_FACTOR = 1.000001;

constexpr double VELOCITY_CONSTRAINTS_EVOLUTION = 2e-2;
constexpr double DEFAULT_TIME_STEP = 2.0e-1;

// Demo split of the DXF outline used by main(): a Cuboid spanning the left half
// of the normalised [-1, 1] design space produced by Loader::alignPoints().
constexpr double DEFAULT_TAG_REGION_MIN = -1.1;
constexpr double DEFAULT_TAG_REGION_SPLIT = 0.0;
constexpr double DEFAULT_TAG_REGION_MAX = 1.1;

// Simulation loop defaults
constexpr unsigned int DEFAULT_NUM_THREADS = 19;
constexpr double DEFAULT_SIMULATION_END_TIME = 14.0;
constexpr unsigned int DEFAULT_SIM_STEPS_PER_FRAME = 5;
constexpr unsigned int DEFAULT_TARGET_FRAMERATE = 60;
constexpr unsigned int DEFAULT_ANTIALIASING_LEVEL = 16;

// Recording defaults
constexpr unsigned int DEFAULT_RECORDING_FRAME_RATE = 20;
constexpr const char* DEFAULT_RECORDING_FRAMES_DIR = "recorded_structure_frames";
constexpr const char* DEFAULT_RECORDING_FRAME_PREFIX = "frame_";
constexpr const char* DEFAULT_ENERGY_LOG_FILENAME = "energy_log.csv";
constexpr const char* DEFAULT_ENERGY_SNAPSHOT_FILENAME = "energy_log_snapshot.png";
constexpr const char* DEFAULT_STRUCTURE_VIDEO_FILENAME = "structure_evolution.mp4";
constexpr const char* ENERGY_PLOT_SCRIPT_FILENAME = "plot_energy_log.py";
constexpr const char* DEFAULT_PYTHON_INTERPRETER = "/bin/python3";

using SpMat = Eigen::SparseMatrix<double>;
using VecXd = Eigen::VectorXd;
using Mat22 = Eigen::Matrix2d;
using Vec2  = Eigen::Vector2d;
