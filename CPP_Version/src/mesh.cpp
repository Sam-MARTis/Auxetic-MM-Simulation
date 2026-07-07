#include "mesh.hpp"

#include "constants.hpp"
#include "utils.hpp"

void Mesh::insertNodesAtPositions(const std::vector<Vec2D> nodePos)
{
    const unsigned int lastIndexOfExistingNodes = nodes.nodesPos.size() - 1;
    const unsigned int numAddedNodes = nodePos.size();
    nodes.numNodes += numAddedNodes;
    nodes.nodesPos.insert(nodes.nodesPos.end(), nodePos.begin(), nodePos.end());
    nodes.constraints.numFree += numAddedNodes;
    std::vector<unsigned int> newFreeIndices;
    newFreeIndices.reserve(numAddedNodes);
    for (unsigned int i = 0; i < numAddedNodes; i++) newFreeIndices.push_back(i+1+lastIndexOfExistingNodes);

    nodes.constraints.freeIndices.insert(nodes.constraints.freeIndices.end(), newFreeIndices.begin(), newFreeIndices.end());

}


void Mesh::initializeBoundaryNodes(const BoundaryMeshBuildInformation& boundaryBuildInfo)
{
    unsigned int indexStart = nodes.nodesPos.size();
    const unsigned int numBoundaryPoints = boundaryBuildInfo.boundary.numBoundaryPoints;
    const std::vector<Vec2D>& boundaryPoints = boundaryBuildInfo.boundary.boundaryPointsPos;
    const bool shouldIncludeCorner = boundaryBuildInfo.pointAtCorners;
    std::vector<Vec2D> pointsToInsert;

    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        std::vector<Vec2D> pointsToInsertThisEdge;
        const Vec2D corner1 = boundaryPoints.at(i);
        const Vec2D corner2 = boundaryPoints.at((i+1)%numBoundaryPoints);
        const Vec2D direction = corner2 - corner1;
        double distance = direction.norm();
        if (!shouldIncludeCorner) distance -= 2*DEFAULT_DISTANCE_FROM_CORNER;
        const unsigned int numNodesToInsert = static_cast<unsigned int>(std::ceil(distance * boundaryBuildInfo.targetDensity))+1;
        pointsToInsertThisEdge.reserve(numNodesToInsert);

        const double distanceBetweenNodes = distance/(numNodesToInsert-1);
        const Vec2D directionUnit = direction.normalized();
        Vec2D startingPoint = corner1;
        if (!shouldIncludeCorner) startingPoint += directionUnit*DEFAULT_DISTANCE_FROM_CORNER;
        for (unsigned int j = 0; j < numNodesToInsert; j++)
        {
            const Vec2D insertPosition = startingPoint + (directionUnit * j*distanceBetweenNodes);
            pointsToInsertThisEdge.push_back(insertPosition);
        }
        BoundaryEdgeInformation thisBoundaryEdgeInfo;
        thisBoundaryEdgeInfo.id = i;
        thisBoundaryEdgeInfo.boundsIndices = std::make_pair(indexStart, indexStart+ numNodesToInsert);
        boundaryEdgeInfo.push_back(thisBoundaryEdgeInfo);
        indexStart += numNodesToInsert;

        pointsToInsert.insert(pointsToInsert.end(), pointsToInsertThisEdge.begin(), pointsToInsertThisEdge.end());


    }
    insertNodesAtPositions(pointsToInsert);

}

void Mesh::initializeInternalEdgesNodes(const InternalEdgesMeshBuildInformation internalEdgeBuildInfo)
{
    const unsigned int numInternalEdges = internalEdgeBuildInfo.internalEdges.numInternalEdges;
    const std::vector<std::pair<Vec2D, Vec2D>>& edgePoints = internalEdgeBuildInfo.internalEdges.internalEdgesPointsPos;
    const bool shouldIncludeCorner = internalEdgeBuildInfo.pointAtCorners;
    std::vector<Vec2D> pointsToInsert;
    for (unsigned int i = 0; i < numInternalEdges; i++)
    {
        std::vector<Vec2D> pointsToInsertThisEdge;
        const std::pair<Vec2D, Vec2D>& edgeCurrent = edgePoints.at(i);
        const Vec2D corner1 = std::get<0>(edgeCurrent);
        const Vec2D corner2 = std::get<1>(edgeCurrent);
        const Vec2D direction = corner2 - corner1;
        const Vec2D normal = Vec2D(-direction.y, direction.x).normalized();
        double distance = direction.norm();
        if (!shouldIncludeCorner) distance -= 2*DEFAULT_DISTANCE_FROM_CORNER;
        const unsigned int numNodesToInsert = static_cast<unsigned int>(std::ceil(distance * internalEdgeBuildInfo.targetDensity))+1;
        pointsToInsertThisEdge.reserve(2*numNodesToInsert);

        const double distanceBetweenNodes = distance/(numNodesToInsert-1);
        const Vec2D directionUnit = direction.normalized();
        Vec2D startingPoint = corner1;
        if (!shouldIncludeCorner) startingPoint += directionUnit*DEFAULT_DISTANCE_FROM_CORNER;
        const Vec2D insertOffsetFromEdge = normal * DEFAULT_INTERNAL_EDGE_INSERT_OFFSET;
        for (unsigned int j = 0; j < numNodesToInsert; j++)
        {
            const Vec2D insertPosition = startingPoint + (directionUnit * j*distanceBetweenNodes);
            pointsToInsertThisEdge.push_back(insertPosition + insertOffsetFromEdge);
            pointsToInsertThisEdge.push_back(insertPosition - insertOffsetFromEdge);
        }

        pointsToInsert.insert(pointsToInsert.end(), pointsToInsertThisEdge.begin(), pointsToInsertThisEdge.end());


    }
    insertNodesAtPositions(pointsToInsert);
}
