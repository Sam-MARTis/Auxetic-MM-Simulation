#pragma once

#include <string>
#include <Eigen/Sparse>
#include <Eigen/Dense>
// #define NDEBUG


constexpr double EPSILON = 1.0e-8;
constexpr double DEFAULT_DISTANCE_FROM_CORNER = 1.0e-4;
constexpr double DEFAULT_INTERNAL_EDGE_INSERT_OFFSET = 1.0e-2;
constexpr double DEFAULT_PARTICLE_RENDER_RADIUS = 1e-2;
constexpr double DEFAULT_SCREEN_WIDTH = 800;
constexpr double DEFAULT_SCREEN_HEIGHT = 800;
constexpr double DEFAULT_SCREEN_SCALE = DEFAULT_SCREEN_WIDTH/2.5;
constexpr unsigned int DEFAULT_NUMBER_INSERT_TRIES = 100;
constexpr unsigned int DEFAULT_NUMBER_INTERNAL_NODES = 500;
constexpr double TARGET_BOUNDARY_NODE_DENSITY = 20.0;
constexpr double TARGET_INTERNAL_EDGE_NODE_DENSITY = 20.0;
constexpr unsigned int DEFAULT_NUMBER_AREA_SAMPLES = 10000;
constexpr unsigned int DEFAULT_NUMBER_NODES_PER_TREE_LEAF = 10;
constexpr unsigned int DEFAULT_NUM_NEIGHBOURS_FOR_EDGES = 55;
constexpr double DEFAULT_DISTANCE_THRESHOLD_FOR_EDGES = 9e-2;
constexpr double DEFAULT_EDGE_RENDER_THICKNESS = 4e-3;
constexpr double DEFAULT_EDGE_STIFFNESS = 1.0;

static constexpr int    MAX_SOLVER_INTERNAL_ITER = 50;
static constexpr double TOL      = 1e-10;
static constexpr double EPSILON2  = 1e-12;

constexpr double DEFAULT_INTERNAL_EDGE_EXPANSIVE_FACTOR = 1.000001;

constexpr double VELOCITY_CONSTRAINTS_EVOLUTION = 1.0e-1;
constexpr double DEFAULT_TIME_STEP = 1.0e-1;

using SpMat = Eigen::SparseMatrix<double>;
using Vec   = Eigen::VectorXd;
using Mat22 = Eigen::Matrix2d;
using Vec2  = Eigen::Vector2d;
