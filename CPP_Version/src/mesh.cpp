#include "mesh.hpp"

#include "constants.hpp"
#include "utils.hpp"

void Mesh::insertNodesAtPositions(const std::vector<Vec2D> nodePos)
{
    const unsigned int lastIndexOfExistingNodes = nodes.nodesPos.size() - 1;
    const unsigned int numAddedNodes = nodePos.size();
    nodes.numNodes += numAddedNodes;
    nodes.nodesPos.insert(nodes.nodesPos.begin(), nodePos.begin(), nodePos.end());
    nodes.constraints.numFree += numAddedNodes;
    std::vector<unsigned int> newFreeIndices;
    newFreeIndices.reserve(numAddedNodes);
    for (unsigned int i = 0; i < numAddedNodes; i++) newFreeIndices.push_back(i+1+lastIndexOfExistingNodes);

    nodes.constraints.freeIndices.insert(nodes.constraints.freeIndices.begin(), newFreeIndices.begin(), newFreeIndices.end());

}


void Mesh::initializeBoundaryNodes(const BoundaryMeshBuildInformation& buildInfo)
{
    const unsigned int numBoundaryPoints = buildInfo.boundary.numBoundaryPoints;
    const std::vector<Vec2D>& boundaryPoints = buildInfo.boundary.boundaryPoints;
    const bool shouldIncludeCorner = buildInfo.pointAtCorners;
    std::vector<Vec2D> pointsToInsert;
    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        std::vector<Vec2D> pointsToInsertThisEdge;
        const Vec2D corner1 = boundaryPoints.at(i);
        const Vec2D corner2 = boundaryPoints.at((i+1)%numBoundaryPoints);
        const Vec2D direction = corner2 - corner1;
        double distance = direction.norm();
        if (!shouldIncludeCorner) distance -= 2*DISTANCE_FROM_CORNER;
        const unsigned int numNodesToInsert = static_cast<unsigned int>(std::ceil(distance * buildInfo.targetDensity))+1;
        pointsToInsertThisEdge.reserve(numNodesToInsert);

        const double distanceBetweenNodes = distance/(numNodesToInsert-1);
        const Vec2D directionUnit = direction.normalized();
        Vec2D startingPoint = corner1;
        if (!shouldIncludeCorner) startingPoint += directionUnit*DISTANCE_FROM_CORNER;
        for (unsigned int j = 0; j < numNodesToInsert; j++)
        {
            const Vec2D insertPosition = startingPoint + (directionUnit * j*distanceBetweenNodes);
            pointsToInsertThisEdge.push_back(insertPosition);
        }

        pointsToInsert.insert(pointsToInsert.end(), pointsToInsertThisEdge.begin(), pointsToInsertThisEdge.end());


    }
    insertNodesAtPositions(pointsToInsert);

}