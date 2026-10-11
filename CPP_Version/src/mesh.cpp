#include "mesh.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_set>

#include "constants.hpp"
#include "utils.hpp"
#include <nanoflann.hpp>
#include <math.h>

namespace
{
    // > 0 means the polygon winds counter clockwise.
    double signedAreaTwice(const std::vector<Vec2D>& polygon)
    {
        double twiceArea = 0.0;
        const std::size_t numPoints = polygon.size();
        for (std::size_t i = 0; i < numPoints; i++)
        {
            const Vec2D& p = polygon.at(i);
            const Vec2D& q = polygon.at((i + 1) % numPoints);
            twiceArea += p[0] * q[1] - q[0] * p[1];
        }
        return twiceArea;
    }

    // Outward normal of a polygon edge: (dy,-dx) when CCW, (-dy,dx) when CW.
    Vec2D outwardNormalOfChord(const Vec2D& direction, const bool counterClockwise)
    {
        const double magnitude = direction.norm();
        if (magnitude <= EPSILON)
        {
            return Vec2D{0.0, 0.0};
        }
        const Vec2D candidate = counterClockwise
            ? Vec2D{direction[1], -direction[0]}
            : Vec2D{-direction[1], direction[0]};
        return candidate / magnitude;
    }
}

void Mesh::resizeNodeAttributes()
{
    const Index numNodes = static_cast<Index>(nodes.nodesPos.size());
    nodes.numNodes = numNodes;
    nodes.edgeIndices.resize(numNodes);
    // resize() preserves existing entries, so tags/normals survive appending.
    nodes.boundaryId.resize(numNodes, NO_BOUNDARY);
    nodes.normals.resize(numNodes, Vec2D{0.0, 0.0});
}

void Mesh::insertNodesAtPositions(const std::vector<Vec2D> nodePos)
{
    if (nodePos.empty()) return;

    nodes.nodesPos.insert(nodes.nodesPos.end(), nodePos.begin(), nodePos.end());
    resizeNodeAttributes();
}

Index Mesh::createBoundaryTag()
{
    const Index boundaryId = numBoundaries++;
    boundaryIdentifiers.emplace_back(boundaryId, std::vector<Index>{});
    boundaryConditions.emplace_back(BoundaryCondition::free());
    return boundaryId;
}

void Mesh::rebuildBoundaryIdentifiers()
{
    for (auto& identifier : boundaryIdentifiers)
    {
        identifier.second.clear();
    }
    numBoundaryNodes = 0;
    for (Index nodeIdx = 0; nodeIdx < nodes.numNodes; nodeIdx++)
    {
        const Index boundaryId = nodes.boundaryId.at(nodeIdx);
        if (boundaryId == NO_BOUNDARY) continue;
        if (boundaryId >= numBoundaries)
        {
            std::cerr << "Mesh::rebuildBoundaryIdentifiers: node " << nodeIdx
                      << " carries unknown boundary id " << boundaryId << std::endl;
            nodes.boundaryId.at(nodeIdx) = NO_BOUNDARY;
            continue;
        }
        boundaryIdentifiers.at(boundaryId).second.push_back(nodeIdx);
        numBoundaryNodes++;
    }
}

Index Mesh::tagNodesByDomain(const Domain& region, Index boundaryId, const bool includeInteriorNodes)
{
    if (boundaryId == NO_BOUNDARY || boundaryId >= numBoundaries)
    {
        boundaryId = createBoundaryTag();
    }

    Index numRetagged = 0;
    for (Index nodeIdx = 0; nodeIdx < nodes.numNodes; nodeIdx++)
    {
        const Index currentTag = nodes.boundaryId.at(nodeIdx);
        if (currentTag == boundaryId) continue;
        // By default only already-tagged boundary nodes are regrouped.
        if (currentTag == NO_BOUNDARY && !includeInteriorNodes) continue;
        if (!region.isInside(nodes.nodesPos.at(nodeIdx))) continue;
        nodes.boundaryId.at(nodeIdx) = boundaryId;
        numRetagged++;
    }

    rebuildBoundaryBookkeeping();
    return boundaryId;
}

void Mesh::pruneBoundarySegments()
{
    std::vector<BoundarySegment> surviving;
    surviving.reserve(boundarySegments.size());
    for (const BoundarySegment& segment : boundarySegments)
    {
        if (segment.node1 >= nodes.numNodes || segment.node2 >= nodes.numNodes) continue;
        if (segment.node1 == segment.node2) continue;
        // Retagging a sub-chain detaches it from the old boundary.
        if (nodes.boundaryId.at(segment.node1) != segment.boundaryId) continue;
        if (nodes.boundaryId.at(segment.node2) != segment.boundaryId) continue;
        surviving.push_back(segment);
    }
    boundarySegments.swap(surviving);
}

void Mesh::inferBoundarySegments(const Index boundaryId)
{
    const std::vector<Index>& taggedNodes = boundaryIdentifiers.at(boundaryId).second;
    if (taggedNodes.size() < 2 || edges.numEdges == 0) return;

    // Candidate chords: spring edges with both endpoints on this tag.
    std::vector<BoundarySegment> candidates;
    for (Index edgeIdx = 0; edgeIdx < edges.numEdges; edgeIdx++)
    {
        const EdgesConnection& connection = edges.edgesConnection.at(edgeIdx);
        if (nodes.boundaryId.at(connection.node1) != boundaryId) continue;
        if (nodes.boundaryId.at(connection.node2) != boundaryId) continue;
        candidates.push_back(BoundarySegment{connection.node1, connection.node2, boundaryId, Vec2D{0.0, 0.0}});
    }
    if (candidates.empty()) return;

    // Keep the two shortest chords per node: a boundary chain is one dimensional,
    // so this drops the edge builder's shortcuts between non-adjacent nodes.
    const auto chordLength = [this](const BoundarySegment& segment)
    {
        return (nodes.nodesPos.at(segment.node1) - nodes.nodesPos.at(segment.node2)).norm();
    };
    std::sort(candidates.begin(), candidates.end(),
              [&chordLength](const BoundarySegment& lhs, const BoundarySegment& rhs)
              {
                  return chordLength(lhs) < chordLength(rhs);
              });

    std::vector<Index> incidentCount(nodes.numNodes, 0);
    std::vector<BoundarySegment> accepted;
    for (const BoundarySegment& candidate : candidates)
    {
        if (chordLength(candidate) <= EPSILON) continue;
        if (incidentCount.at(candidate.node1) >= 2) continue;
        if (incidentCount.at(candidate.node2) >= 2) continue;
        incidentCount.at(candidate.node1)++;
        incidentCount.at(candidate.node2)++;
        accepted.push_back(candidate);
    }

    // Orient away from the tag centroid; ambiguous for internal boundaries, where
    // the caller should supply BoundaryCondition::normalOverride instead.
    Vec2D centroid{0.0, 0.0};
    for (const Index nodeIdx : taggedNodes)
    {
        centroid += nodes.nodesPos.at(nodeIdx);
    }
    centroid /= static_cast<double>(taggedNodes.size());

    for (BoundarySegment& segment : accepted)
    {
        const Vec2D chord = nodes.nodesPos.at(segment.node2) - nodes.nodesPos.at(segment.node1);
        Vec2D normal = outwardNormalOfChord(chord, false);
        const Vec2D midpoint = (nodes.nodesPos.at(segment.node1) + nodes.nodesPos.at(segment.node2)) * 0.5;
        if (Vec2D::dot(normal, midpoint - centroid) < 0.0)
        {
            normal *= -1.0;
        }
        segment.outwardNormal = normal;
        boundarySegments.push_back(segment);
    }
}

void Mesh::rebuildNormalsForBoundary(const Index boundaryId, const Vec2D& normalOverride)
{
    const std::vector<Index>& taggedNodes = boundaryIdentifiers.at(boundaryId).second;
    if (taggedNodes.empty()) return;

    for (const Index nodeIdx : taggedNodes)
    {
        nodes.normals.at(nodeIdx) = Vec2D{0.0, 0.0};
    }

    const double overrideMagnitude = normalOverride.norm();
    if (overrideMagnitude > EPSILON)
    {
        const Vec2D unitNormal = normalOverride / overrideMagnitude;
        for (const Index nodeIdx : taggedNodes)
        {
            nodes.normals.at(nodeIdx) = unitNormal;
        }
        return;
    }

    for (const BoundarySegment& segment : boundarySegments)
    {
        if (segment.boundaryId != boundaryId) continue;
        nodes.normals.at(segment.node1) += segment.outwardNormal;
        nodes.normals.at(segment.node2) += segment.outwardNormal;
    }

    Vec2D centroid{0.0, 0.0};
    for (const Index nodeIdx : taggedNodes)
    {
        centroid += nodes.nodesPos.at(nodeIdx);
    }
    centroid /= static_cast<double>(taggedNodes.size());

    for (const Index nodeIdx : taggedNodes)
    {
        const Vec2D accumulated = nodes.normals.at(nodeIdx);
        const double magnitude = accumulated.norm();
        if (magnitude > EPSILON)
        {
            nodes.normals.at(nodeIdx) = accumulated / magnitude;
            continue;
        }
        // Lone tagged node: fall back to the radial direction.
        const Vec2D radial = nodes.nodesPos.at(nodeIdx) - centroid;
        const double radialMagnitude = radial.norm();
        nodes.normals.at(nodeIdx) =
            radialMagnitude > EPSILON ? (radial / radialMagnitude) : Vec2D{0.0, 0.0};
    }
}

void Mesh::computeBoundaryNodeNormals()
{
    std::fill(nodes.normals.begin(), nodes.normals.end(), Vec2D{0.0, 0.0});
    for (Index boundaryId = 0; boundaryId < numBoundaries; boundaryId++)
    {
        rebuildNormalsForBoundary(boundaryId, boundaryConditions.at(boundaryId).normalOverride);
    }
}

void Mesh::buildNodeEdgeAdjacency()
{
    for (Index nodeIdx = 0; nodeIdx < nodes.numNodes; nodeIdx++)
    {
        nodes.edgeIndices.at(nodeIdx).clear();
    }
    for (Index edgeIdx = 0; edgeIdx < edges.numEdges; edgeIdx++)
    {
        const EdgesConnection& connection = edges.edgesConnection.at(edgeIdx);
        if (connection.node1 < nodes.numNodes) nodes.edgeIndices.at(connection.node1).push_back(edgeIdx);
        if (connection.node2 < nodes.numNodes) nodes.edgeIndices.at(connection.node2).push_back(edgeIdx);
    }
}

void Mesh::classifyBoundaryEdges()
{
    boundaryEdgesIndices.clear();
    domainEdgesIndices.clear();
    for (Index edgeIdx = 0; edgeIdx < edges.numEdges; edgeIdx++)
    {
        EdgeProperties& property = edges.edgesProperty.at(edgeIdx);
        const EdgesConnection& connection = edges.edgesConnection.at(edgeIdx);
        const Index tag1 = nodes.boundaryId.at(connection.node1);
        const Index tag2 = nodes.boundaryId.at(connection.node2);
        if (tag1 != NO_BOUNDARY && tag1 == tag2)
        {
            property.boundaryId = tag1;
            boundaryEdgesIndices.push_back(edgeIdx);
        }
        else
        {
            property.boundaryId = NO_BOUNDARY;
            domainEdgesIndices.push_back(edgeIdx);
        }
    }
}

void Mesh::rebuildBoundaryBookkeeping()
{
    rebuildBoundaryIdentifiers();
    if (edges.numEdges > 0) buildNodeEdgeAdjacency();
    pruneBoundarySegments();
    for (Index boundaryId = 0; boundaryId < numBoundaries; boundaryId++)
    {
        const bool hasSegments = std::any_of(
            boundarySegments.begin(), boundarySegments.end(),
            [boundaryId](const BoundarySegment& segment) { return segment.boundaryId == boundaryId; });
        if (!hasSegments) inferBoundarySegments(boundaryId);
    }
    computeBoundaryNodeNormals();
    if (edges.numEdges > 0) classifyBoundaryEdges();
    bookkeepingDirty = false;
}

void Mesh::ensureBoundaryBookkeeping()
{
    if (bookkeepingDirty) rebuildBoundaryBookkeeping();
}

void Mesh::setBoundaryCondition(const Index boundaryId, const BoundaryCondition& condition)
{
    if (boundaryId >= numBoundaries)
    {
        std::cerr << "Mesh::setBoundaryCondition: unknown boundary id " << boundaryId
                  << " (numBoundaries=" << numBoundaries << ")" << std::endl;
        return;
    }
    boundaryConditions.at(boundaryId) = condition;
    if (bookkeepingDirty)
    {
        rebuildBoundaryBookkeeping();
    }
    else
    {
        rebuildNormalsForBoundary(boundaryId, condition.normalOverride);
    }
}

const BoundaryCondition& Mesh::boundaryCondition(const Index boundaryId) const
{
    return boundaryConditions.at(boundaryId);
}

const std::vector<Index>& Mesh::boundaryNodes(const Index boundaryId) const
{
    return boundaryIdentifiers.at(boundaryId).second;
}

bool Mesh::isBoundaryNode(const Index nodeIndex) const
{
    return nodes.boundaryId.at(nodeIndex) != NO_BOUNDARY;
}

void Mesh::logBoundarySummary() const
{
    std::cout << "Boundaries: " << numBoundaries
              << ", tagged nodes: " << numBoundaryNodes
              << ", boundary edges: " << boundaryEdgesIndices.size()
              << ", domain edges: " << domainEdgesIndices.size() << std::endl;
    for (Index boundaryId = 0; boundaryId < numBoundaries; boundaryId++)
    {
        const BoundaryCondition& condition = boundaryConditions.at(boundaryId);
        std::cout << "  boundary " << boundaryId
                  << " | nodes " << boundaryIdentifiers.at(boundaryId).second.size()
                  << " | " << boundaryConditionTypeName(condition.type)
                  << " | velocity " << condition.velocity;
        if (!boundaryIdentifiers.at(boundaryId).second.empty())
        {
            const Index probe = boundaryIdentifiers.at(boundaryId).second.front();
            const Vec2D normal = nodes.normals.at(probe);
            std::cout << " | normal@" << probe << " (" << normal[0] << ", " << normal[1] << ")";
        }
        std::cout << std::endl;
    }
}

// ---------------------------------------------------------------------------
// Mesh generation
// ---------------------------------------------------------------------------

Index Mesh::initializeBoundaryNodes(const BoundaryMeshBuildInformation& boundaryBuildInfo)
{
    const Index thisBoundaryNumPoints = static_cast<Index>(boundaryBuildInfo.boundary.numPoints());
    const std::vector<Vec2D>& boundaryPoints = boundaryBuildInfo.boundary.vertices();
    const bool shouldIncludeCorner = boundaryBuildInfo.pointAtCorners;
    const Index boundaryId = createBoundaryTag();

    const bool counterClockwise = thisBoundaryNumPoints >= 3 && signedAreaTwice(boundaryPoints) > 0.0;

    std::vector<Vec2D> pointsToInsert;
    std::vector<Index> nodesPerEdge(thisBoundaryNumPoints, 0);
    std::vector<Vec2D> normalPerEdge(thisBoundaryNumPoints, Vec2D{0.0, 0.0});

    for (Index i = 0; i < thisBoundaryNumPoints; i++)
    {
        std::vector<Vec2D> pointsToInsertThisEdge;
        const Vec2D corner1 = boundaryPoints.at(i);
        const Vec2D corner2 = boundaryPoints.at((i + 1) % thisBoundaryNumPoints);

        const Vec2D direction = corner2 - corner1;
        double distance = direction.norm();
        if (!shouldIncludeCorner) distance -= 2 * DEFAULT_DISTANCE_FROM_CORNER;
        if (distance <= EPSILON) continue;

        const Index numNodesToInsert =
            static_cast<Index>(std::ceil(distance * boundaryBuildInfo.targetDensity)) + 1;
        if (numNodesToInsert < 2) continue;

        pointsToInsertThisEdge.reserve(numNodesToInsert);
        const double distanceBetweenNodes = distance / (numNodesToInsert - 1);
        const Vec2D directionUnit = direction.normalized();
        Vec2D startingPoint = corner1;
        if (!shouldIncludeCorner) startingPoint += directionUnit * DEFAULT_DISTANCE_FROM_CORNER;

        for (Index j = 0; j < numNodesToInsert; j++)
        {
            pointsToInsertThisEdge.push_back(startingPoint + directionUnit * (j * distanceBetweenNodes));
        }

        nodesPerEdge.at(i) = numNodesToInsert;
        normalPerEdge.at(i) = outwardNormalOfChord(direction, counterClockwise);
        pointsToInsert.insert(pointsToInsert.end(), pointsToInsertThisEdge.begin(), pointsToInsertThisEdge.end());
    }

    const Index firstNewIndex = static_cast<Index>(nodes.nodesPos.size());
    insertNodesAtPositions(pointsToInsert);

    Index cursor = firstNewIndex;
    for (Index i = 0; i < thisBoundaryNumPoints; i++)
    {
        const Index count = nodesPerEdge.at(i);
        if (count == 0) continue;
        const Vec2D edgeNormal = normalPerEdge.at(i);
        for (Index k = 0; k < count; k++)
        {
            nodes.boundaryId.at(cursor + k) = boundaryId;
        }
        for (Index k = 0; k + 1 < count; k++)
        {
            boundarySegments.push_back(BoundarySegment{cursor + k, cursor + k + 1, boundaryId, edgeNormal});
        }
        cursor += count;
    }

    // Owned copy: the caller's BoundaryMeshBuildInformation may be a temporary.
    ownedDomainFilters.push_back(std::make_unique<PointHull>(boundaryBuildInfo.boundary));
    domainFilters.push_back(ownedDomainFilters.back().get());
    domainFilterParity.push_back(true);
    bookkeepingDirty = true;
    return boundaryId;
}

void Mesh::initializeInternalEdgesNodes(const InternalEdgesMeshBuildInformation internalEdgeBuildInfo)
{
    if (domainFilters.empty())
    {
        std::cerr << "Boundary nodes must be initialized before internal edge nodes" << std::endl;
        exit(1);
    }

    const Index numInternalEdges = internalEdgeBuildInfo.internalEdges.numInternalEdges;
    const std::vector<std::pair<Vec2D, Vec2D>>& edgePoints = internalEdgeBuildInfo.internalEdges.internalEdgesPointsPos;
    const bool shouldIncludeCorner = internalEdgeBuildInfo.pointAtCorners;
    std::vector<Vec2D> pointsToInsert;
    for (Index i = 0; i < numInternalEdges; i++)
    {
        std::vector<Vec2D> pointsToInsertThisEdge;
        const std::pair<Vec2D, Vec2D>& edgeCurrent = edgePoints.at(i);
        segmentDiscontinuities.push_back(edgeCurrent);
        numSegmentDiscontinuities++;

        const Vec2D corner1 = std::get<0>(edgeCurrent);
        const Vec2D corner2 = std::get<1>(edgeCurrent);
        const Vec2D direction = corner2 - corner1;
        const Vec2D normal = Vec2D{-direction[1], direction[0]}.normalized();
        double distance = direction.norm();
        if (!shouldIncludeCorner) distance -= 2 * DEFAULT_DISTANCE_FROM_CORNER;
        if (distance <= EPSILON) continue;
        const Index numNodesToInsert =
            static_cast<Index>(std::ceil(distance * internalEdgeBuildInfo.targetDensity)) + 1;
        if (numNodesToInsert < 2) continue;
        pointsToInsertThisEdge.reserve(2 * numNodesToInsert);

        const double distanceBetweenNodes = distance / (numNodesToInsert - 1);
        const Vec2D directionUnit = direction.normalized();
        Vec2D startingPoint = corner1;
        if (!shouldIncludeCorner) startingPoint += directionUnit * DEFAULT_DISTANCE_FROM_CORNER;
        const Vec2D insertOffsetFromEdge = normal * DEFAULT_INTERNAL_EDGE_INSERT_OFFSET;
        for (Index j = 0; j < numNodesToInsert; j++)
        {
            const Vec2D insertPosition = startingPoint + directionUnit * (j * distanceBetweenNodes);
            const Vec2D p1 = insertPosition + insertOffsetFromEdge;
            const Vec2D p2 = insertPosition - insertOffsetFromEdge;
            const Vec2D points[2] = {p1, p2};
            for (const Vec2D& p : points)
            {
                bool isValidPoint = true;
                for (Index filterIdx = 0; filterIdx < static_cast<Index>(domainFilters.size()); filterIdx++)
                {
                    if (domainFilters.at(filterIdx)->isInside(p) != domainFilterParity.at(filterIdx))
                    {
                        isValidPoint = false;
                        break;
                    }
                }
                if (isValidPoint) pointsToInsertThisEdge.push_back(p);
            }
        }

        pointsToInsert.insert(pointsToInsert.end(), pointsToInsertThisEdge.begin(), pointsToInsertThisEdge.end());
    }
    insertNodesAtPositions(pointsToInsert);
}

void Mesh::addInternalNodes_Random(const Index numNodesToAdd)
{
    std::vector<Vec2D> pointsToAdd;
    for (Index i = 0; i < numNodesToAdd; i++)
    {
        Vec2D randPos;
        for (Index iter = 0; iter < 1; iter++)
        {
            randPos[0] = randValueMine(-1.0, 1.0);
            randPos[1] = randValueMine(-1.0, 1.0);
            bool pointSatisfied = true;
            for (Index filterIdx = 0; filterIdx < static_cast<Index>(domainFilters.size()); filterIdx++)
            {
                if (domainFilters.at(filterIdx)->isInside(randPos) != domainFilterParity.at(filterIdx))
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

double Mesh::calculateAreaOfMeshObject(const Index numSamples)
{
    const double areaOfPhysicsSpace = 4.0;
    Index numInside = 0;
    for (Index i = 0; i < numSamples; i++)
    {
        bool isValid = true;
        const Vec2D randPos = Vec2D{randValueMine(-1.0, 1.0), randValueMine(-1.0, 1.0)};
        for (Index filterIdx = 0; filterIdx < static_cast<Index>(domainFilters.size()); filterIdx++)
        {
            if (domainFilters.at(filterIdx)->isInside(randPos) != domainFilterParity.at(filterIdx))
            {
                isValid = false;
                break;
            }
        }
        if (isValid) numInside++;
    }
    areaOfDxf = areaOfPhysicsSpace * (static_cast<double>(numInside)) / (static_cast<double>(numSamples));
    areaCalculationDone = true;
    return areaOfDxf;
}

void Mesh::addInternalNodes_Uniform(const int numNodesToTargetAdd)
{
    if (!areaCalculationDone)
    {
        std::cerr << "Area calculation not done. Uniform addition of nodes not allowed" << std::endl;
        exit(1);
    }
    std::vector<Vec2D> pointsToInsert;
    const Index numNodesToAdd = static_cast<Index>((4.0 * numNodesToTargetAdd) / areaOfDxf);
    const Index numNodesWillAddPerDim = static_cast<Index>(sqrt(static_cast<double>(numNodesToAdd)));
    const double ds = 2.0 / (numNodesWillAddPerDim + 1.0);
    const double x_start = ds - 1.0;
    const double y_start = ds - 1.0;
    for (Index i = 0; i < numNodesWillAddPerDim; i++)
    {
        const double x = x_start + i * ds;
        for (Index j = 0; j < numNodesWillAddPerDim; j++)
        {
            const double y = y_start + j * ds;
            const Vec2D pointTest = Vec2D{x, y};

            bool isValidPoint = true;
            for (Index filterIdx = 0; filterIdx < static_cast<Index>(domainFilters.size()); filterIdx++)
            {
                if (domainFilters.at(filterIdx)->isInside(pointTest) != domainFilterParity.at(filterIdx))
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

void Mesh::findNearestNeighbours(const KDTree& tree, const Vec2D& queryPoint, const Index K,
                                 std::vector<size_t>& indices, std::vector<double>& distances)
{
    indices.resize(K);
    distances.resize(K);
    nanoflann::KNNResultSet<double> result(K);
    result.init(indices.data(), distances.data());
    tree.findNeighbors(result, queryPoint.data.data());
}

void Mesh::createEdges(const Index K, const double distanceTheshold)
{
    PointCloud cloud{.pointsPos = nodes.nodesPos};
    KDTree tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(DEFAULT_NUMBER_NODES_PER_TREE_LEAF));
    std::vector<size_t> indices(K);
    std::vector<double> distances(K);
    const double distanceThresholdSq = distanceTheshold * distanceTheshold;
    std::unordered_set<EdgeKey, EdgeHash> edgeSet;
    assert(nodes.numNodes == nodes.nodesPos.size());
    for (Index node1Idx = 0; node1Idx < nodes.numNodes; node1Idx++)
    {
        const Vec2D& queryPoint = nodes.nodesPos.at(node1Idx);
        findNearestNeighbours(tree, queryPoint, K, indices, distances);
        for (Index edgeIdx = 0; edgeIdx < K; edgeIdx++)
        {
            const Index node2Idx = static_cast<Index>(indices.at(edgeIdx));
            const double distanceSq = distances.at(edgeIdx);
            // Discard the KD-tree self-hit (zero length, infinite stiffness).
            if (node2Idx == node1Idx) continue;
            if (distanceSq > distanceThresholdSq) continue;
            if (doesEdgeIntersectDiscontinuity(node1Idx, node2Idx)) continue;
            EdgeKey edge(node1Idx, node2Idx);
            if (!(edgeSet.insert(edge).second)) continue;
            const Vec2D& p1 = nodes.nodesPos.at(edge.a);
            const Vec2D& p2 = nodes.nodesPos.at(edge.b);
            const double naturalLength = (p1 - p2).norm();
            if (naturalLength <= EPSILON) continue;
            edges.edgesConnection.push_back(EdgesConnection{edge.a, edge.b});
            edges.edgesProperty.push_back(EdgeProperties{
                .natural_length = naturalLength,
                .stiffness = DEFAULT_EDGE_STIFFNESS / naturalLength,
                .cached_extension = naturalLength,
            });
            // Derived geometry caches.
            edges.vertices.push_back(std::array<Vec2D, 2>{p1, p2});
            edges.edgeEvaluationPoints.push_back((p1 + p2) * 0.5);
            edges.numEdges++;
        }
    }
    rebuildBoundaryBookkeeping();
}

bool Mesh::doesEdgeIntersectDiscontinuity(const Index node1Idx, const Index node2Idx)
{
    for (Index segmentIdx = 0; segmentIdx < numSegmentDiscontinuities; segmentIdx++)
    {
        const std::pair<Vec2D, Vec2D>& segment = segmentDiscontinuities.at(segmentIdx);
        if (VecMath::segmentsIntersects2D(nodes.nodesPos.at(node1Idx), nodes.nodesPos.at(node2Idx),
                                          std::get<0>(segment), std::get<1>(segment)))
        {
            return true;
        }
    }
    return false;
}

void Mesh::createFlatVec()
{
    nodePosFlat = VecXd(nodes.numNodes * 2);
    for (Index nodeIdx = 0; nodeIdx < nodes.numNodes; nodeIdx++)
    {
        nodePosFlat[2 * nodeIdx] = nodes.nodesPos.at(nodeIdx)[0];
        nodePosFlat[2 * nodeIdx + 1] = nodes.nodesPos.at(nodeIdx)[1];
    }
}

void Mesh::cacheEdgesExtensionFromFlatVec()
{
    for (Index i = 0; i < edges.numEdges; i++)
    {
        EdgeProperties& property = edges.edgesProperty.at(i);
        const EdgesConnection& connection = edges.edgesConnection.at(i);
        const double x1 = nodePosFlat[2 * connection.node1];
        const double y1 = nodePosFlat[2 * connection.node1 + 1];
        const double x2 = nodePosFlat[2 * connection.node2];
        const double y2 = nodePosFlat[2 * connection.node2 + 1];
        const double dx = x1 - x2;
        const double dy = y1 - y2;
        property.cached_extension = sqrt(dx * dx + dy * dy);

        const Vec2D p1{x1, y1};
        const Vec2D p2{x2, y2};
        edges.vertices.at(i) = std::array<Vec2D, 2>{p1, p2};
        edges.edgeEvaluationPoints.at(i) = (p1 + p2) * 0.5;
    }
}

void Mesh::syncVecWithFlatVecData()
{
    for (Index i = 0; i < nodes.numNodes; i++)
    {
        nodes.nodesPos.at(i)[0] = nodePosFlat[2 * i];
        nodes.nodesPos.at(i)[1] = nodePosFlat[2 * i + 1];
    }
}
