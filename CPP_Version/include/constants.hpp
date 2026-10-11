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
constexpr double DEFAULT_BOUNDARY_NORMAL_RENDER_LENGTH = 4e-2;

static constexpr int    MAX_SOLVER_INTERNAL_ITER = 5;
static constexpr double TOL      = 1e-2;
static constexpr double EPSILON2  = 1e-12;

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
