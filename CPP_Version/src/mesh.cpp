#include "mesh.hpp"

#include <iostream>

#include "constants.hpp"
#include "utils.hpp"
#include <nanoflann.hpp>

void Mesh::insertNodesAtPositions(const std::vector<Vec2D> nodePos)
{
    const unsigned int numAddedNodes = nodePos.size();
    if (numAddedNodes == 0) return;

    const unsigned int firstNewIndex = static_cast<unsigned int>(nodes.nodesPos.size());
    nodes.numNodes += numAddedNodes;
    nodes.nodesPos.insert(nodes.nodesPos.end(), nodePos.begin(), nodePos.end());
    nodes.constraints.numFree += numAddedNodes;
    std::vector<unsigned int> newFreeIndices;
    newFreeIndices.reserve(numAddedNodes);
    for (unsigned int i = 0; i < numAddedNodes; i++) newFreeIndices.push_back(firstNewIndex + i);

    nodes.constraints.freeIndices.insert(nodes.constraints.freeIndices.end(), newFreeIndices.begin(), newFreeIndices.end());
    nodes.numNodes = static_cast<unsigned int>(nodes.nodesPos.size());

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
        segmentDiscontinuities.push_back(std::make_pair(corner1, corner2));
        numSegmentDiscontinuities ++;

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
        thisBoundaryEdgeInfo.id = highestBoundaryID;
        highestBoundaryID++;
        thisBoundaryEdgeInfo.boundsIndices = std::make_pair(indexStart, indexStart+ numNodesToInsert);
        boundaryEdgeInfo.push_back(thisBoundaryEdgeInfo);
        indexStart += numNodesToInsert;

        pointsToInsert.insert(pointsToInsert.end(), pointsToInsertThisEdge.begin(), pointsToInsertThisEdge.end());


    }
    insertNodesAtPositions(pointsToInsert);
    boundaries.push_back(boundaryBuildInfo.boundary);
    numBoundaries++;
    boundaryIncludeExclude.push_back(true);
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
        segmentDiscontinuities.push_back(edgeCurrent);
        numSegmentDiscontinuities++;

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

void Mesh::addInternalNodes_Random(const unsigned int numNodesToAdd)
{
    std::vector<Vec2D> pointsToAdd;
    for (unsigned int i = 0; i < numNodesToAdd; i++)
    {
        Vec2D randPos;
        for (unsigned int iter = 0; iter< 1; iter++)
        {
            randPos.x = randValueMine(-1.0, 1.0);
            randPos.y = randValueMine(-1.0, 1.0);
            bool pointSatisfied = true;
            for (unsigned int boundaryIdx = 0; boundaryIdx < numBoundaries; boundaryIdx++)
            {
                if (boundaries.at(boundaryIdx).isInsideBoundary(randPos) != boundaryIncludeExclude.at(boundaryIdx))
                {
                    pointSatisfied = false;

                    break;
                }
            }
            if (pointSatisfied)
            {
                pointsToAdd.push_back(randPos);
                break;
            }
        }

    }
    insertNodesAtPositions(pointsToAdd);

}

double Mesh::calculateAreaOfMeshObject(const unsigned int numSamples)
{
    const double areaOfPhysicsSpace = 4.0;
    unsigned int numInside = 0;
    for (unsigned int i = 0; i < numSamples; i++)
    {
        bool isValid = true;
        for (unsigned int boundaryIdx = 0; boundaryIdx < numBoundaries; boundaryIdx++)
        {
            const Vec2D randPos = Vec2D(randValueMine(-1.0, 1.0), randValueMine(-1.0, 1.0));
            if (boundaries.at(boundaryIdx).isInsideBoundary(Vec2D(randPos)) != boundaryIncludeExclude.at(boundaryIdx))
            {
                isValid = false;
                break;
            }
        }
        if (isValid) numInside++;
    }
    areaOfDxf = areaOfPhysicsSpace * (static_cast<double>(numInside))/(static_cast<double>(numSamples));
    areaCalculationDone = true;
    return areaOfDxf;
}

void Mesh::addInternalNodes_Uniform(const int numNodesToTargetAdd)
{
    boundaries.at(0).printBoundaryPoints();
    if (!areaCalculationDone)
    {
        std::cerr << "Area calculation not done. Uniform addition of nodes not allowed"<<std::endl;
        exit(1);
    }
    std::vector<Vec2D> pointsToInsert;
    const unsigned int numNodesToAdd = static_cast<unsigned int>((4.0 * numNodesToTargetAdd)/areaOfDxf);
    const unsigned int numNodesWillAddPerDim = static_cast<unsigned int>(sqrt(static_cast<double>(numNodesToAdd)));
    const double ds = 2.0/(numNodesWillAddPerDim+1.0);
    double x_start = ds-1.0;
    double y_start = ds-1.0;
    std::cout<<"Num boundaries "<<numBoundaries<<std::endl;
    for (unsigned int i = 0; i<numNodesWillAddPerDim; i++)
    {
        double x = x_start + i*ds;
        for (unsigned int j = 0; j<numNodesWillAddPerDim; j++)
        {
            double y = y_start + j*ds;
            const Vec2D pointTest = Vec2D(x,y);

            bool isValidPoint = true;
            for (unsigned int boundaryIdx = 0; boundaryIdx < numBoundaries; boundaryIdx++)
            {
                if (boundaries.at(boundaryIdx).isInsideBoundary(pointTest) != boundaryIncludeExclude.at(boundaryIdx))
                {
                    isValidPoint = false;
                    break;
                }
            }
            if (isValidPoint)
            {
                pointsToInsert.push_back(pointTest);
            }
        }
    }
    insertNodesAtPositions(pointsToInsert);
}

void Mesh::findNearestNeighbours(const KDTree& tree, const Vec2D& queryPoint, unsigned int K, std::vector<size_t>& indices, std::vector<double>& distances)
{
    indices.resize(K);
    distances.resize(K);
    nanoflann::KNNResultSet<double> result(K);
    result.init(indices.data(), distances.data());
    tree.findNeighbors(result, &queryPoint.x);
}

void Mesh::createEdges(const unsigned int K, const double distanceTheshold)
{
    PointCloud cloud{
    .pointsPos = nodes.nodesPos};
    KDTree tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(DEFAULT_NUMBER_NODES_PER_TREE_LEAF));
    std::vector<size_t> indices(K);
    std::vector<double> distances(K);
    const double distanceThresholdSq = distanceTheshold*distanceTheshold;
    std::unordered_set<EdgeKey, EdgeHash> edgeSet;
    assert(nodes.numNodes == nodes.nodesPos.size());
    for (unsigned int node1Idx=0; node1Idx<nodes.numNodes; node1Idx++)
    {
        const Vec2D& queryPoint = nodes.nodesPos.at(node1Idx);
        findNearestNeighbours(tree, queryPoint, K, indices, distances);
        for (unsigned int edgeIdx = 0; edgeIdx<K; edgeIdx++)
        {
            const unsigned int node2Idx = indices.at(edgeIdx);
            const double distanceSq = distances.at(edgeIdx);
            if (distanceSq > distanceThresholdSq) continue;
            EdgeKey edge(node1Idx, node2Idx);
            if (!(edgeSet.insert(edge).second)) continue;
            edges.edgesConnection.push_back({edge.a, edge.b});
            edges.edgesProperty.push_back({
                .natural_length = (nodes.nodesPos.at(edge.a)  - nodes.nodesPos.at(edge.b)).norm(),
                .stiffness = DEFAULT_EDGE_STIFFNESS,
                .cached_extension = 0.0
            });
            edges.numEdges++;
        }
    }
}

bool Mesh::doesEdgeIntersectDiscontinuity(const EdgesConnection& edgeToTest)
{
    for (unsigned int segmentIdx = 0; segmentIdx < numSegmentDiscontinuities; segmentIdx++)
    {
        const std::pair<Vec2D, Vec2D>& segment = segmentDiscontinuities.at(segmentIdx);
        if (segmentsIntersects2D(nodes.nodesPos.at(edgeToTest.node1), nodes.nodesPos.at(edgeToTest.node2), std::get<0>(segment), std::get<1>(segment)))
        {
            return true;
        }
    }
    return false;
}
