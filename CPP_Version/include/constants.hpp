#pragma once
#include <string>
// #define NDEBUG

constexpr double EPSILON = 1.0e-8;
constexpr double DEFAULT_DISTANCE_FROM_CORNER = 1.0e-4;
constexpr double DEFAULT_INTERNAL_EDGE_INSERT_OFFSET = 1.0e-3;
constexpr double DEFAULT_PARTICLE_RENDER_RADIUS = 1e-2;
constexpr double DEFAULT_EDGE_RENDER_THICKNESS = 2e-2;
constexpr double DEFAULT_SCREEN_WIDTH = 800;
constexpr double DEFAULT_SCREEN_HEIGHT = 800;
constexpr double DEFAULT_SCREEN_SCALE = DEFAULT_SCREEN_WIDTH/2.5;
constexpr unsigned int DEFAULT_NUMBER_INSERT_TRIES = 100;
constexpr unsigned int DEFAULT_NUMBER_INTERNAL_NODES = 1000;
constexpr unsigned int DEFAULT_NUMBER_AREA_SAMPLES = 10000;
constexpr unsigned int DEFAULT_NUMBER_NODES_PER_TREE_LEAF = 10;
constexpr double DEFAULT_EDGE_STIFFNESS = 1.0;