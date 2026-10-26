#pragma once
#include <cstddef>
#include <vector>
#include <cmath>

#include <tuple>
#include "utils.hpp"
#include "domains.hpp"
#include <string>

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
    unsigned int numBoundaryLines;
    unsigned int numCutPairs;
public:
    // Closed, ordered polygon assembled from the DXF's BOUNDARY lines.
    std::vector<Vec2D> boundaryNodesPos;
    std::vector<std::pair<Vec2D, Vec2D>> internalEdgesNodesPos;
    Loader(const std::string& dxfFilepath);
    void alignPoints();
    [[nodiscard]] std::size_t numBoundaryPolygonPoints() const { return boundaryNodesPos.size(); }
};
