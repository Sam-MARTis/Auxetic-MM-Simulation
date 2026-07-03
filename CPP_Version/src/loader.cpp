#include "loader.hpp"
#include "utils.hpp"

Boundary::Boundary(const std::vector<Vec2D>& points)
{
    this->boundaryPoints.reserve(points.size());
    // Add the first point repeat in the end
    this->boundaryPoints.insert(this->boundaryPoints.end(), points.begin(), points.end());

}


