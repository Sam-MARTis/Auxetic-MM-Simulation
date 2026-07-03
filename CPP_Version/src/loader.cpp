#include "loader.hpp"

#include "constants.hpp"
#include "utils.hpp"

Boundary::Boundary(const std::vector<Vec2D>& points)
{
    // Copy points and set the count used by isInsideBoundary
    this->boundaryPoints = points;
    this->numBoundaryPoints = static_cast<unsigned int>(points.size());
}

bool Boundary::isInsideBoundary(const Vec2D& testPoint) const
{
    const Vec2D rayDirection = Vec2D(1.0, 0.0);
    unsigned int numIntersections = 0;
    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        const Vec2D& boundaryPoint = boundaryPoints[i];
        const Vec2D segmentDirection = boundaryPoints[(i+1)%numBoundaryPoints] - boundaryPoints[i];
        // const double denom = rayDirection.x*segmentDirection.y - rayDirection.y*segmentDirection.x;
        const double denom = crossVec2D(rayDirection, segmentDirection);
        if (std::abs(denom) <= EPSILON) continue;
        const double t = crossVec2D(boundaryPoint - testPoint, segmentDirection)/denom;
        if (t < 0.0) continue;
        const double u = crossVec2D(boundaryPoint - testPoint, rayDirection)/denom;
        if (u>=0.0 && u < 1.0) numIntersections++;

    }
    if (numIntersections%2 == 0) return false;
    return true;
}
