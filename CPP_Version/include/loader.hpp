#pragma once
#include <vector>
#include <cmath>

#include <tuple>
#include "utils.hpp"
#include <string>

// struct Segment
// {
//     Point2D startPos;
//     Point2D endPos;
//     Segment(Point2D startPos, Point2D endPos){this->startPos=startPos; this->endPos=endPos;};
// };

class Boundary // Might be worth considering if I should replace this with a Domain class
{
public:
    Boundary(const std::vector<Vec2D>& points);
    unsigned int numBoundaryPoints;
    std::vector<Vec2D> boundaryPointsPos; // Probably better to store indices, eh?
    // Wait is this just the defining boundary points? Oh shit, I remember now.
    [[nodiscard]] bool isInsideBoundary(const Vec2D& testPoint) const;
    void printBoundaryPoints() const;
    void testPoint(const Vec2D& testPoint) const; // Gives a description of the point. Stuff liek whether it is inside or outside.
};

struct BoundaryMeshBuildInformation
{
    double targetDensity=10.0;
    Boundary boundary;
    bool pointAtCorners = false;
};

struct InternalEdgesDescription
{

    InternalEdgesDescription(const std::vector<std::pair<Vec2D, Vec2D>>& points);
    unsigned int numInternalEdges;
    std::vector<std::pair<Vec2D, Vec2D>> internalEdgesPointsPos;
};

struct InternalEdgesMeshBuildInformation
{
    double targetDensity=10.0;
    InternalEdgesDescription internalEdges;
    bool pointAtCorners = false;
};

class Loader
{
private:
    unsigned int numBoundaryPoints;
    unsigned int numCutPairs;
public:
    std::vector<Vec2D> boundaryNodesPos;
    std::vector<std::pair<Vec2D, Vec2D>> internalEdgesNodesPos;
    Loader(const std::string& dxfFilepath);
    void alignPoints();
};

