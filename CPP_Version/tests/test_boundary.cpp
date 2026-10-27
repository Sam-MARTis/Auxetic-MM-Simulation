// Host-only tests for the mesh data structure, boundary tagging and boundary
// conditions. Headless; exit code 0 means all checks passed.

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "boundary.hpp"
#include "constants.hpp"
#include "domains.hpp"
#include "mesh.hpp"
#include "solid.hpp"

namespace
{
    int numChecks = 0;
    int numFailures = 0;

    void expect(const bool condition, const std::string& label)
    {
        numChecks++;
        if (condition)
        {
            std::cout << "  ok   " << label << std::endl;
        }
        else
        {
            numFailures++;
            std::cout << "  FAIL " << label << std::endl;
        }
    }

    void expectNear(const double actual, const double expected, const double tolerance,
                    const std::string& label)
    {
        const bool ok = std::abs(actual - expected) <= tolerance;
        numChecks++;
        if (ok)
        {
            std::cout << "  ok   " << label << std::endl;
        }
        else
        {
            numFailures++;
            std::cout << "  FAIL " << label << " (actual " << actual
                      << ", expected " << expected << ", tol " << tolerance << ")" << std::endl;
        }
    }

    PointHull makeSquare(const double halfSize)
    {
        return PointHull(std::vector<Vec2D>{
            {-halfSize, -halfSize}, {halfSize, -halfSize},
            {halfSize, halfSize}, {-halfSize, halfSize}});
    }

    // A tagged square sheet with a uniform interior grid and springs.
    void buildSquareSheet(Mesh& mesh, const Index numInteriorNodes)
    {
        BoundaryMeshBuildInformation boundaryInfo{
            .targetDensity = TARGET_BOUNDARY_NODE_DENSITY,
            .boundary = makeSquare(1.0),
        };
        mesh.initializeBoundaryNodes(boundaryInfo);
        mesh.calculateAreaOfMeshObject(DEFAULT_NUMBER_AREA_SAMPLES);
        mesh.addInternalNodes_Uniform(static_cast<int>(numInteriorNodes));
        mesh.createEdges(DEFAULT_NUM_NEIGHBOURS_FOR_EDGES, DEFAULT_DISTANCE_THRESHOLD_FOR_EDGES);
    }

    void perturbInteriorNodes(Mesh& mesh, const double amplitude)
    {
        for (Index nodeIdx = 0; nodeIdx < mesh.nodes.numNodes; nodeIdx++)
        {
            if (mesh.nodes.boundaryId.at(nodeIdx) != NO_BOUNDARY) continue;
            const double x = mesh.nodes.nodesPos.at(nodeIdx)[0];
            const double y = mesh.nodes.nodesPos.at(nodeIdx)[1];
            mesh.nodes.nodesPos.at(nodeIdx) = Vec2D{x + amplitude * y, y - amplitude * x};
        }
    }

    Vec2D positionOf(const Mesh& mesh, const Index nodeIdx)
    {
        return Vec2D{mesh.nodePosFlat[2 * nodeIdx], mesh.nodePosFlat[2 * nodeIdx + 1]};
    }

    void testPointHull()
    {
        std::cout << "[PointHull::isInside]" << std::endl;
        const PointHull square = makeSquare(1.0);
        expect(square.isInside(Vec2D{0.5, 0.5}), "interior point is inside");
        expect(!square.isInside(Vec2D{1.5, 0.5}), "exterior point is outside");
        expect(square.isInside(Vec2D{0.0, 0.0}), "centre is inside");

        const PointHull triangle(std::vector<Vec2D>{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}});
        expect(triangle.isInside(Vec2D{0.2, 0.2}), "triangle interior point is inside");
        expect(!triangle.isInside(Vec2D{0.8, 0.8}), "triangle exterior point is outside");
    }

    void testAutoTaggingAndNormals()
    {
        std::cout << "[auto tagging, adjacency and outward normals]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);

        expect(mesh.numBoundaries == 1, "one boundary tag is created per polygon");
        expect(mesh.boundaryConditions.size() == 1, "every tag has a condition slot");
        expect(mesh.boundaryConditions.at(0).type == BoundaryConditionType::FREE,
               "tags default to FREE");
        expect(mesh.numBoundaryNodes > 0, "boundary nodes were tagged");
        expect(mesh.boundaryIdentifiers.at(0).first == 0, "boundaryIdentifiers is indexed by id");
        expect(mesh.boundaryIdentifiers.at(0).second.size() == mesh.numBoundaryNodes,
               "boundaryIdentifiers[0] lists every tagged node");

        Index numTagged = 0;
        double maxUnitNormalError = 0.0;
        double minOutwardDot = std::numeric_limits<double>::max();
        for (Index nodeIdx = 0; nodeIdx < mesh.nodes.numNodes; nodeIdx++)
        {
            if (mesh.nodes.boundaryId.at(nodeIdx) == NO_BOUNDARY) continue;
            numTagged++;
            const Vec2D normal = mesh.nodes.normals.at(nodeIdx);
            maxUnitNormalError = std::max(maxUnitNormalError, std::abs(normal.norm() - 1.0));
            minOutwardDot = std::min(minOutwardDot,
                                     Vec2D::dot(normal, mesh.nodes.nodesPos.at(nodeIdx)));
        }
        expect(numTagged == mesh.numBoundaryNodes, "per-node tags agree with the id map");
        expect(maxUnitNormalError < 1e-9, "every tagged node carries a unit normal");
        expect(minOutwardDot > 0.9, "every normal points away from the square's centre");

        expect(mesh.boundaryEdgesIndices.size() + mesh.domainEdgesIndices.size() == mesh.edges.numEdges,
               "every edge is classified as boundary or domain");
        bool boundaryEdgesValid = true;
        for (const Index edgeIdx : mesh.boundaryEdgesIndices)
        {
            const EdgesConnection& connection = mesh.edges.edgesConnection.at(edgeIdx);
            if (mesh.nodes.boundaryId.at(connection.node1) != mesh.nodes.boundaryId.at(connection.node2))
                boundaryEdgesValid = false;
            if (mesh.edges.edgesProperty.at(edgeIdx).boundaryId == NO_BOUNDARY)
                boundaryEdgesValid = false;
        }
        expect(boundaryEdgesValid, "boundary edges join two nodes of one tag and carry that tag");

        std::size_t adjacencyEntries = 0;
        for (const std::vector<Index>& incidentEdges : mesh.nodes.edgeIndices)
            adjacencyEntries += incidentEdges.size();
        expect(adjacencyEntries == 2 * mesh.edges.numEdges,
               "adjacency lists exactly two endpoints per edge");
        bool adjacencyValid = true;
        for (Index edgeIdx = 0; edgeIdx < mesh.edges.numEdges; edgeIdx++)
        {
            const EdgesConnection& connection = mesh.edges.edgesConnection.at(edgeIdx);
            const std::vector<Index>& first = mesh.nodes.edgeIndices.at(connection.node1);
            const std::vector<Index>& second = mesh.nodes.edgeIndices.at(connection.node2);
            if (std::find(first.begin(), first.end(), edgeIdx) == first.end()) adjacencyValid = false;
            if (std::find(second.begin(), second.end(), edgeIdx) == second.end()) adjacencyValid = false;
        }
        expect(adjacencyValid, "every edge appears in both endpoint adjacency lists");
    }

    void testDerivedEdgeGeometryCaches()
    {
        std::cout << "[FVM-parity derived edge geometry caches]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        mesh.createFlatVec();
        mesh.cacheEdgesExtensionFromFlatVec();

        double maxMidpointError = 0.0;
        double maxVertexError = 0.0;
        double maxLengthError = 0.0;
        for (Index edgeIdx = 0; edgeIdx < mesh.edges.numEdges; edgeIdx++)
        {
            const EdgesConnection& connection = mesh.edges.edgesConnection.at(edgeIdx);
            const Vec2D p1 = mesh.nodes.nodesPos.at(connection.node1);
            const Vec2D p2 = mesh.nodes.nodesPos.at(connection.node2);
            const Vec2D expectedMidpoint = (p1 + p2) * 0.5;
            maxMidpointError = std::max(maxMidpointError,
                (mesh.edges.edgeEvaluationPoints.at(edgeIdx) - expectedMidpoint).norm());
            maxVertexError = std::max(maxVertexError,
                (mesh.edges.vertices.at(edgeIdx)[0] - p1).norm());
            maxLengthError = std::max(maxLengthError,
                std::abs(mesh.edges.edgesProperty.at(edgeIdx).cached_extension - (p1 - p2).norm()));
        }
        expect(maxVertexError < 1e-12, "vertices mirror the endpoint positions");
        expect(maxMidpointError < 1e-12, "edgeEvaluationPoints hold the midpoints");
        expect(maxLengthError < 1e-12, "cached_extension holds the current length");
    }

    void testDomainRetagging()
    {
        std::cout << "[Domain-based retagging (FVM identifyBoundaryEdgesAndCells analogue)]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        const Index originalCount = mesh.numBoundaryNodes;

        const Cuboid leftHalf(Vec2D{-1.1, -1.1}, Vec2D{0.0, 1.1});
        const Index newId = mesh.tagNodesByDomain(leftHalf);

        expect(newId == 1, "a fresh tag id is handed out");
        expect(mesh.numBoundaries == 2, "numBoundaries grew to two");
        expect(!mesh.boundaryNodes(1).empty(), "the region captured boundary nodes");
        expect(mesh.boundaryNodes(0).size() + mesh.boundaryNodes(1).size() == originalCount,
               "the two tags partition the original boundary nodes");
        expect(mesh.numBoundaryNodes == originalCount, "retagging preserves the tagged node count");

        bool partitionValid = true;
        for (Index nodeIdx = 0; nodeIdx < mesh.nodes.numNodes; nodeIdx++)
        {
            const Index tag = mesh.nodes.boundaryId.at(nodeIdx);
            if (tag == NO_BOUNDARY) continue;
            const bool inside = leftHalf.isInside(mesh.nodes.nodesPos.at(nodeIdx));
            if ((tag == newId) != inside) partitionValid = false;
        }
        expect(partitionValid, "tag 1 is exactly the set of boundary nodes inside the region");

        double minOutwardDot = std::numeric_limits<double>::max();
        double maxUnitNormalError = 0.0;
        Index worstNode = NO_BOUNDARY;
        for (const Index nodeIdx : mesh.boundaryNodes(newId))
        {
            const Vec2D normal = mesh.nodes.normals.at(nodeIdx);
            maxUnitNormalError = std::max(maxUnitNormalError, std::abs(normal.norm() - 1.0));
            const double outwardDot = Vec2D::dot(normal, mesh.nodes.nodesPos.at(nodeIdx));
            if (outwardDot < minOutwardDot)
            {
                minOutwardDot = outwardDot;
                worstNode = nodeIdx;
            }
        }
        expect(maxUnitNormalError < 1e-9, "inferred domain-tag normals are unit vectors");
        if (minOutwardDot <= 0.9)
        {
            const Vec2D p = mesh.nodes.nodesPos.at(worstNode);
            const Vec2D n = mesh.nodes.normals.at(worstNode);
            std::cout << "       worst normal node " << worstNode << " at (" << p[0] << ", " << p[1]
                      << ") n=(" << n[0] << ", " << n[1] << ") dot=" << minOutwardDot << std::endl;
            for (const BoundarySegment& segment : mesh.boundarySegments)
            {
                if (segment.node1 == worstNode || segment.node2 == worstNode)
                {
                    std::cout << "         segment b" << segment.boundaryId
                              << " (" << segment.node1 << "," << segment.node2 << ") n=("
                              << segment.outwardNormal[0] << ", " << segment.outwardNormal[1] << ")"
                              << std::endl;
                }
            }
            const auto incident = mesh.nodes.edgeIndices.at(worstNode);
            std::cout << "         incident springs:";
            for (const Index e : incident)
            {
                const EdgesConnection& c = mesh.edges.edgesConnection.at(e);
                const Index other = (c.node1 == worstNode) ? c.node2 : c.node1;
                std::cout << " " << other << "(tag " << mesh.nodes.boundaryId.at(other) << ")";
            }
            std::cout << std::endl;
        }
        expect(minOutwardDot > 0.9, "inferred domain-tag normals point outward");

        // Retagging into an existing id reuses it instead of creating a new one.
        const Cuboid rightOfSplit(Vec2D{0.05, -1.1}, Vec2D{1.1, 1.1});
        const Index reused = mesh.tagNodesByDomain(rightOfSplit, 0);
        expect(reused == 0, "passing an existing id reuses that tag");
        expect(mesh.numBoundaries == 2, "reusing an id does not create a new tag");
        bool reuseValid = true;
        for (const Index nodeIdx : mesh.boundaryNodes(0))
        {
            if (mesh.nodes.nodesPos.at(nodeIdx)[0] <= 0.0) reuseValid = false;
        }
        expect(reuseValid, "tag 0 still holds only nodes right of the split");
        expect(mesh.boundaryNodes(0).size() + mesh.boundaryNodes(1).size() == originalCount,
               "the partition is still complete after reuse");
    }

    void testInteriorNodePromotion()
    {
        std::cout << "[promoting interior nodes into a tag]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        const Index boundaryNodesBefore = mesh.numBoundaryNodes;

        const Cuboid centre(Vec2D{-0.3, -0.3}, Vec2D{0.3, 0.3});
        const Index patchId = mesh.tagNodesByDomain(centre, NO_BOUNDARY, true);

        expect(mesh.numBoundaryNodes > boundaryNodesBefore,
               "interior nodes were promoted into the new tag");
        expect(!mesh.boundaryNodes(patchId).empty(), "the internal patch is non-empty");
        bool allInside = true;
        for (const Index nodeIdx : mesh.boundaryNodes(patchId))
        {
            if (nodeIdx >= mesh.nodes.numNodes) allInside = false;
            if (!centre.isInside(mesh.nodes.nodesPos.at(nodeIdx))) allInside = false;
        }
        expect(allInside, "every promoted node really lies inside the region");
    }

    void testFixedPinsBoundary()
    {
        std::cout << "[FIXED: nodes pinned, interior relaxes]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        // A small shear keeps the perturbed state in the well-behaved regime.
        perturbInteriorNodes(mesh, 0.01);
        mesh.setBoundaryCondition(0, BoundaryCondition::fixed());
        mesh.createFlatVec();
        Solid solid(mesh);

        const std::vector<Index> tagged = mesh.boundaryNodes(0);
        expect(solid.numFixedNodesCount() == tagged.size(),
               "FIXED eliminates exactly the tagged nodes' DOFs");
        expect(solid.numConstrainedNodes() == 0, "FIXED adds no plane-constraint rows");

        mesh.cacheEdgesExtensionFromFlatVec();
        const double energyBefore = solid.computeElasticEnergyFromExtensionCache();
        expect(energyBefore > 0.0, "the perturbed state carries some elastic energy");

        for (int step = 0; step < 5; step++) solid.timeStep(DEFAULT_TIME_STEP);

        double maxBoundaryDrift = 0.0;
        for (const Index nodeIdx : tagged)
        {
            maxBoundaryDrift = std::max(maxBoundaryDrift,
                (positionOf(mesh, nodeIdx) - mesh.nodes.nodesPos.at(nodeIdx)).norm());
        }
        expect(maxBoundaryDrift < 1e-12, "boundary nodes never move");

        double maxInteriorMove = 0.0;
        for (Index nodeIdx = 0; nodeIdx < mesh.nodes.numNodes; nodeIdx++)
        {
            if (mesh.nodes.boundaryId.at(nodeIdx) != NO_BOUNDARY) continue;
            maxInteriorMove = std::max(maxInteriorMove,
                (positionOf(mesh, nodeIdx) - mesh.nodes.nodesPos.at(nodeIdx)).norm());
        }
        expect(maxInteriorMove > 1e-6, "interior nodes relaxed away from the perturbed state");
        expect(solid.computeElasticEnergyFromExtensionCache() < energyBefore,
               "relaxation lowered the elastic energy");
    }

    void testFixedNormalHoldsThePlane()
    {
        std::cout << "[FIXED_NORMAL: plane held, tangential sliding still allowed]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        // The right half drives outwards and drags the left half along its planes.
        const Index pinnedHalf = mesh.tagNodesByDomain(Cuboid(Vec2D{-1.1, -1.1}, Vec2D{0.0, 1.1}));
        mesh.setBoundaryCondition(0, BoundaryCondition::displaceNormal(0.05));
        mesh.setBoundaryCondition(pinnedHalf, BoundaryCondition::fixedNormal());
        mesh.createFlatVec();
        Solid solid(mesh);

        const std::vector<Index> tagged = mesh.boundaryNodes(pinnedHalf);
        std::vector<Vec2D> start;
        std::vector<Vec2D> normals;
        start.reserve(tagged.size());
        normals.reserve(tagged.size());
        for (const Index nodeIdx : tagged)
        {
            start.push_back(mesh.nodes.nodesPos.at(nodeIdx));
            normals.push_back(mesh.nodes.normals.at(nodeIdx));
        }

        expect(solid.numFixedNodesCount() == 0, "FIXED_NORMAL does not eliminate DOFs");
        Index expectedPlaneRows = 0;
        for (Index nodeIdx = 0; nodeIdx < mesh.nodes.numNodes; nodeIdx++)
            expectedPlaneRows += static_cast<Index>(mesh.nodes.constraintNormals.at(nodeIdx).size());
        expect(solid.numConstrainedNodes() == expectedPlaneRows,
               "one plane row per incident boundary plane");
        expect(solid.numConstrainedNodes() >= mesh.numBoundaryNodes,
               "every tagged node has at least one plane row");

        for (int step = 0; step < 5; step++) solid.timeStep(DEFAULT_TIME_STEP);

        double maxPlaneDrift = 0.0;
        double maxTangentialMove = 0.0;
        double minNormalAdvance = std::numeric_limits<double>::max();
        for (std::size_t k = 0; k < tagged.size(); k++)
        {
            const Vec2D delta = positionOf(mesh, tagged[k]) - start[k];
            const double normalAdvance = Vec2D::dot(normals[k], delta);
            maxPlaneDrift = std::max(maxPlaneDrift, std::abs(normalAdvance));
            maxTangentialMove = std::max(maxTangentialMove,
                                         (delta - normals[k] * normalAdvance).norm());
        }
        for (const Index nodeIdx : mesh.boundaryNodes(0))
        {
            minNormalAdvance = std::min(minNormalAdvance,
                Vec2D::dot(mesh.nodes.normals.at(nodeIdx),
                           positionOf(mesh, nodeIdx) - mesh.nodes.nodesPos.at(nodeIdx)));
        }
        expect(maxPlaneDrift < 1e-8, "the pinned half's normal offset stayed constant");
        expect(maxTangentialMove > 1e-6, "pinned-half nodes still slid along their planes");
        expect(minNormalAdvance > 0.0, "the driven half still advanced along its outward normal");
    }

    void testDisplaceNormalMovesOutward()
    {
        std::cout << "[DISPLACE_NORMAL: move normally outward at a prescribed velocity]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        const double velocity = 0.05;
        mesh.setBoundaryCondition(0, BoundaryCondition::displaceNormal(velocity));
        mesh.createFlatVec();
        Solid solid(mesh);

        const std::vector<Index> tagged = mesh.boundaryNodes(0);
        std::vector<Vec2D> start;
        std::vector<Vec2D> normals;
        start.reserve(tagged.size());
        normals.reserve(tagged.size());
        for (const Index nodeIdx : tagged)
        {
            start.push_back(mesh.nodes.nodesPos.at(nodeIdx));
            normals.push_back(mesh.nodes.normals.at(nodeIdx));
        }

        expect(solid.numFixedNodesCount() == 0, "DISPLACE_NORMAL does not eliminate DOFs");
        expect(solid.numConstrainedNodes() >= tagged.size(),
               "DISPLACE_NORMAL adds at least one plane row per tagged node");

        const int numSteps = 4;
        for (int step = 0; step < numSteps; step++) solid.timeStep(DEFAULT_TIME_STEP);
        const double elapsed = numSteps * DEFAULT_TIME_STEP;

        double minAdvance = std::numeric_limits<double>::max();
        double maxAdvanceError = 0.0;
        for (std::size_t k = 0; k < tagged.size(); k++)
        {
            const double advance = Vec2D::dot(normals[k], positionOf(mesh, tagged[k]) - start[k]);
            minAdvance = std::min(minAdvance, advance);
            maxAdvanceError = std::max(maxAdvanceError, std::abs(advance - velocity * elapsed));
        }
        expect(minAdvance > 0.0, "every tagged node advanced along its outward normal");
        expect(maxAdvanceError < 1e-8, "the normal advance equals velocity * time");
    }

    void testInwardVelocityIsSigned()
    {
        std::cout << "[DISPLACE_NORMAL: negative velocity moves inward]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        mesh.setBoundaryCondition(0, BoundaryCondition::displaceNormal(-0.05));
        mesh.createFlatVec();
        Solid solid(mesh);

        const std::vector<Index> tagged = mesh.boundaryNodes(0);
        std::vector<Vec2D> start;
        for (const Index nodeIdx : tagged) start.push_back(mesh.nodes.nodesPos.at(nodeIdx));

        for (int step = 0; step < 3; step++) solid.timeStep(DEFAULT_TIME_STEP);

        double maxAdvance = -std::numeric_limits<double>::max();
        for (std::size_t k = 0; k < tagged.size(); k++)
        {
            const double advance =
                Vec2D::dot(mesh.nodes.normals.at(tagged[k]), positionOf(mesh, tagged[k]) - start[k]);
            maxAdvance = std::max(maxAdvance, advance);
        }
        expect(maxAdvance < 0.0, "every tagged node retreated along its outward normal");
    }

    void testFreeHasNoConstraints()
    {
        std::cout << "[FREE: no constraints at all]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        mesh.setBoundaryCondition(0, BoundaryCondition::free());
        mesh.createFlatVec();
        Solid solid(mesh);

        // All-FREE has rigid-body modes and a singular Hessian, so it is not
        // solved here; only the structural statement is checked.
        expect(solid.numFixedNodesCount() == 0, "FREE eliminates no DOFs");
        expect(solid.numConstrainedNodes() == 0, "FREE adds no plane rows");
    }

    void testFreeBoundaryMovesWhileFixedDoesNot()
    {
        std::cout << "[FREE vs FIXED on two tags of the same mesh]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        perturbInteriorNodes(mesh, 0.05);

        const Index leftId = mesh.tagNodesByDomain(Cuboid(Vec2D{-1.1, -1.1}, Vec2D{0.0, 1.1}));
        mesh.setBoundaryCondition(0, BoundaryCondition::free());
        mesh.setBoundaryCondition(leftId, BoundaryCondition::fixed());
        mesh.createFlatVec();
        Solid solid(mesh);

        expect(solid.numFixedNodesCount() == mesh.boundaryNodes(leftId).size(),
               "only the pinned tag is eliminated");

        std::vector<Vec2D> leftStart;
        for (const Index nodeIdx : mesh.boundaryNodes(leftId))
            leftStart.push_back(mesh.nodes.nodesPos.at(nodeIdx));
        std::vector<Vec2D> rightStart;
        for (const Index nodeIdx : mesh.boundaryNodes(0))
            rightStart.push_back(mesh.nodes.nodesPos.at(nodeIdx));

        for (int step = 0; step < 5; step++) solid.timeStep(DEFAULT_TIME_STEP);

        double maxPinnedMove = 0.0;
        for (std::size_t k = 0; k < leftStart.size(); k++)
        {
            maxPinnedMove = std::max(maxPinnedMove,
                (positionOf(mesh, mesh.boundaryNodes(leftId)[k]) - leftStart[k]).norm());
        }
        expect(maxPinnedMove < 1e-12, "the pinned tag did not move");

        double maxFreeMove = 0.0;
        for (std::size_t k = 0; k < rightStart.size(); k++)
        {
            maxFreeMove = std::max(maxFreeMove,
                (positionOf(mesh, mesh.boundaryNodes(0)[k]) - rightStart[k]).norm());
        }
        expect(maxFreeMove > 1e-6, "the FREE tag moved freely");
    }

    void testConditionCannotBeSilentlyDropped()
    {
        std::cout << "[normal override and unresolved-normal reporting]" << std::endl;
        Mesh mesh;
        buildSquareSheet(mesh, DEFAULT_NUMBER_INTERNAL_NODES);
        const Vec2D override_{0.0, 1.0};
        mesh.setBoundaryCondition(0, BoundaryCondition::fixedNormal(override_));

        bool allOverridden = true;
        for (const Index nodeIdx : mesh.boundaryNodes(0))
        {
            const Vec2D normal = mesh.nodes.normals.at(nodeIdx);
            if ((normal - override_).norm() > 1e-12) allOverridden = false;
        }
        expect(allOverridden, "normalOverride replaces the geometry-derived normal");

        mesh.setBoundaryCondition(0, BoundaryCondition::displaceNormal(0.0, Vec2D{2.0, 0.0}));
        bool normalizedOverride = true;
        for (const Index nodeIdx : mesh.boundaryNodes(0))
        {
            if ((mesh.nodes.normals.at(nodeIdx) - Vec2D{1.0, 0.0}).norm() > 1e-12)
                normalizedOverride = false;
        }
        expect(normalizedOverride, "a non-unit override is normalised");
    }

    void testTorsionalElement()
    {
        std::cout << "[torsional bond-angle element]" << std::endl;

        SpringMeshDescription description;
        description.nodePositions = {
            Vec2D{0.0, 0.0}, Vec2D{1.0, 0.0}, Vec2D{0.35, 0.9}};
        description.springs = {
            EdgesConnection{0, 1}, EdgesConnection{1, 2}, EdgesConnection{2, 0}};
        description.triangles = {std::array<Index, 3>{0, 1, 2}};
        description.numTriangles = 1;

        Mesh mesh;
        mesh.initializeFromSpringMesh(description);

        expect(mesh.edges.numEdges == 3, "one triangle gives three springs");
        expect(mesh.torsions.numTorsions == 3,
               "one triangle gives three torsional corners");
        expect(mesh.numNearCollinearTorsions == 0,
               "the reference triangle is not near-collinear");

        const auto hasSpring = [&mesh](const Index a, const Index b)
        {
            for (Index e = 0; e < mesh.edges.numEdges; e++)
            {
                const EdgesConnection& c = mesh.edges.edgesConnection.at(e);
                if ((c.node1 == a && c.node2 == b) || (c.node1 == b && c.node2 == a))
                    return true;
            }
            return false;
        };

        bool allLegsAreSprings = true;
        bool stiffnessFormula = true;
        for (Index t = 0; t < mesh.torsions.numTorsions; t++)
        {
            const TorsionConnection& c = mesh.torsions.torsionConnection.at(t);
            const TorsionProperties& p = mesh.torsions.torsionProperty.at(t);
            if (!hasSpring(c.nodeI, c.nodeApex) || !hasSpring(c.nodeApex, c.nodeK))
                allLegsAreSprings = false;
            const double kA = mesh.edges.edgesProperty.at(p.legEdgeA).stiffness;
            const double kB = mesh.edges.edgesProperty.at(p.legEdgeB).stiffness;
            const double expected = DEFAULT_TORSIONAL_STIFFNESS
                                  * 0.5 * (kA + kB)
                                  * p.natural_length_a * p.natural_length_b;
            if (std::abs(p.stiffness - expected) > 1e-12) stiffnessFormula = false;
        }
        expect(allLegsAreSprings, "every torsion leg is a spring");
        expect(stiffnessFormula, "torsion stiffness follows the edge-normalised formula");

        mesh.createFlatVec();
        Solid solid(mesh);

        // Perturb away from the reference so the gradient is non-zero.
        for (Index i = 0; i < mesh.nodes.numNodes; i++)
        {
            mesh.nodePosFlat[2 * i]     += 0.04 * static_cast<double>(i + 1);
            mesh.nodePosFlat[2 * i + 1] -= 0.03 * static_cast<double>(i + 1);
        }
        solid.assembleEnergyDerivatives();

        const VecXd g = solid.gradient();
        const SpMat K = solid.hessian();
        const int numDofs = static_cast<int>(g.size());
        expect(numDofs == 6, "three free nodes give six DOFs");

        const VecXd x0 = mesh.nodePosFlat;
        const double h = 1e-5;

        double maxGradientError = 0.0;
        for (int p = 0; p < numDofs; p++)
        {
            mesh.nodePosFlat = x0;
            mesh.nodePosFlat[p] += h;
            const double plus = solid.computeTotalPotentialEnergy();
            mesh.nodePosFlat = x0;
            mesh.nodePosFlat[p] -= h;
            const double minus = solid.computeTotalPotentialEnergy();
            const double fd = (plus - minus) / (2.0 * h);
            maxGradientError = std::max(maxGradientError, std::abs(fd - g[p]));
        }
        mesh.nodePosFlat = x0;
        expect(maxGradientError < 1e-6, "torsion gradient matches finite differences");

        double maxHessianError = 0.0;
        for (int p = 0; p < numDofs; p++)
        {
            for (int q = 0; q < numDofs; q++)
            {
                const auto energyAt = [&](const double dp, const double dq)
                {
                    mesh.nodePosFlat = x0;
                    mesh.nodePosFlat[p] += dp;
                    mesh.nodePosFlat[q] += dq;
                    return solid.computeTotalPotentialEnergy();
                };
                const double fd = (energyAt(h, h) - energyAt(h, -h)
                                   - energyAt(-h, h) + energyAt(-h, -h))
                                  / (4.0 * h * h);
                maxHessianError = std::max(maxHessianError, std::abs(fd - K.coeff(p, q)));
            }
        }
        mesh.nodePosFlat = x0;
        expect(maxHessianError < 1e-4, "torsion Hessian matches finite differences");

        double maxAsymmetry = 0.0;
        for (int p = 0; p < numDofs; p++)
            for (int q = 0; q < numDofs; q++)
                maxAsymmetry = std::max(maxAsymmetry, std::abs(K.coeff(p, q) - K.coeff(q, p)));
        expect(maxAsymmetry < 1e-12, "torsion Hessian is symmetric");

        VecXd translationX = VecXd::Zero(numDofs);
        VecXd translationY = VecXd::Zero(numDofs);
        for (int i = 0; i < numDofs / 2; i++)
        {
            translationX[2 * i] = 1.0;
            translationY[2 * i + 1] = 1.0;
        }
        expect((K * translationX).norm() < 1e-10 && (K * translationY).norm() < 1e-10,
               "torsion Hessian annihilates translations");

        mesh.createFlatVec();
        expect(solid.computeTotalPotentialEnergy() < 1e-12,
               "the reference triangle has zero potential energy");
    }
}

int main()
{
    std::cout << "=== spring mesh / boundary condition tests ===" << std::endl;
    testPointHull();
    testAutoTaggingAndNormals();
    testDerivedEdgeGeometryCaches();
    testDomainRetagging();
    testInteriorNodePromotion();
    testFixedPinsBoundary();
    testFixedNormalHoldsThePlane();
    testDisplaceNormalMovesOutward();
    testInwardVelocityIsSigned();
    testFreeHasNoConstraints();
    testFreeBoundaryMovesWhileFixedDoesNot();
    testConditionCannotBeSilentlyDropped();
    testTorsionalElement();

    std::cout << "=== " << (numChecks - numFailures) << "/" << numChecks
              << " checks passed ===" << std::endl;
    return numFailures == 0 ? 0 : 1;
}
