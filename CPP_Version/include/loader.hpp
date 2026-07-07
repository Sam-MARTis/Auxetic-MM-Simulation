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

class Boundary
{
public:
    Boundary(const std::vector<Vec2D>& points);
    unsigned int numBoundaryPoints;
    std::vector<Vec2D> boundaryPointsPos;
    [[nodiscard]] bool isInsideBoundary(const Vec2D& testPoint) const;
};

struct BoundaryMeshBuildInformation
{
    double targetDensity=10.0;
    Boundary boundary;
    bool pointAtCorners = false;
};

class InternalEdges
{
    public:
    InternalEdges(const std::vector<std::pair<Vec2D, Vec2D>>& points);
    unsigned int numInternalEdges;
    std::vector<std::pair<Vec2D, Vec2D>> internalEdgesPointsPos;

};

struct InternalEdgesMeshBuildInformation
{
    double targetDensity=10.0;
    InternalEdges internalEdges;
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

