#pragma once
#include <vector>
#include <cmath>

#include <tuple>
#include "utils.hpp"
#include "domains.hpp"
#include <string>

// struct Segment
// {
//     Point2D startPos;
//     Point2D endPos;
//     Segment(Point2D startPos, Point2D endPos){this->startPos=startPos; this->endPos=endPos;};
// };

struct BoundaryMeshBuildInformation
{
    double targetDensity=10.0;
    PointHull boundary;
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

