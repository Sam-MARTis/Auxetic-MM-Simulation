#pragma once
#include <vector>
#include <cmath>

#include "utils.hpp"


// struct Segment
// {
//     Point2D startPos;
//     Point2D endPos;
//     Segment(Point2D startPos, Point2D endPos){this->startPos=startPos; this->endPos=endPos;};
// };

class Boundary
{
public:
    Boundary(const std::vector<Vec2D>& points);
    unsigned int numBoundaryPoints;
    std::vector<Vec2D> boundaryPointsPos;
    bool isInsideBoundary(const Vec2D& testPoint) const;
};

struct BoundaryMeshBuildInformation
{
    double targetDensity;
    Boundary boundary;
    bool pointAtCorners = false;
    // BoundaryMeshBuildInformation();

};

    