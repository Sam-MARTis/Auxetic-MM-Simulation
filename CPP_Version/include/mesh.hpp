#pragma once

#include <array>
#include <functional>
#include <memory>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <nanoflann.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "boundary.hpp"
#include "constants.hpp"
#include "domains.hpp"
#include "loader.hpp"
#include "utils.hpp"

struct PointCloud
{
    std::vector<Vec2D>& pointsPos;
    size_t kdtree_get_point_count() const { return pointsPos.size(); };
    double kdtree_get_pt(const size_t idx, int dim) const { return dim==0?pointsPos[idx][0]:pointsPos[idx][1]; };
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false;}
};
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>;

// FVM EdgeProperty::connectingCells analogue.
struct EdgesConnection
{
    Index node1 = NO_BOUNDARY;
    Index node2 = NO_BOUNDARY;
};

// FVM EdgeProperty analogue; parallel to Edges::edgesConnection.
struct EdgeProperties
{
    double natural_length = 0.0;
    double stiffness = 0.0;
    double cached_extension = 0.0;   // current length, refreshed from nodePosFlat
    Index boundaryId = NO_BOUNDARY;  // boundary tag, set only for boundary chords
};

// A ready-made spring network, e.g. derived from a Gmsh mesh.
struct SpringMeshDescription
{
    std::vector<Vec2D> nodePositions;
    std::vector<EdgesConnection> springs;
    std::vector<BoundarySegment> boundarySegments;         // outer boundary chords
    std::vector<std::pair<Vec2D, Vec2D>> cutSegments;
    // Diagnostics filled in by the mesher.
    Index numTriangles = 0;
    Index numBoundaryEdges = 0;
    Index numCutBoundaryEdges = 0;
};

// Cells and edges follow the FVM reference layout: struct-of-structs with
// parallel arrays addressed by index. See DeveloperReference.md section 1.
class Mesh
{
public:
    // One entry per edge in every array, addressed by edge index.
    struct Edges
    {
        Index numEdges = 0;
        std::vector<EdgeProperties> edgesProperty;             // M
        std::vector<EdgesConnection> edgesConnection;          // M
        std::vector<std::array<Vec2D, 2>> vertices;            // M, derived cache
        std::vector<Vec2D> edgeEvaluationPoints;               // M, derived cache
    };

    // A spring network has no cells: the node is the control volume.
    struct Nodes
    {
        Index numNodes = 0;
        std::vector<Vec2D> nodesPos;                           // N
        std::vector<std::vector<Index>> edgeIndices;           // N, incident edges
        std::vector<Index> boundaryId;                         // N, NO_BOUNDARY if interior
        std::vector<Vec2D> normals;                            // N, outward unit normal when tagged
        // N, one entry per distinct incident boundary plane: a straight edge gives
        // one, a corner gives one per edge, so the corner cannot slide.
        std::vector<std::vector<Vec2D>> constraintNormals;
    };

    Nodes nodes;
    Edges edges;
    VecXd nodePosFlat;

    // ---- boundary bookkeeping -------------------------------------------------
    Index numBoundaries = 0;
    Index numBoundaryNodes = 0;
    // Id -> tagged node indices. Ids are dense: boundaryIdentifiers[i].first == i.
    std::vector<std::pair<Index, std::vector<Index>>> boundaryIdentifiers;
    std::vector<Index> boundaryEdgesIndices;
    std::vector<Index> domainEdgesIndices;
    std::vector<BoundarySegment> boundarySegments;
    std::vector<BoundaryCondition> boundaryConditions;         // one per boundary id

    // Material region used while seeding nodes, not boundary tags. Mesh-owned, so
    // passing a temporary BoundaryMeshBuildInformation is safe.
    double areaOfDxf = 4.0;

    // ---- mesh generation ------------------------------------------------------
    void insertNodesAtPositions(const std::vector<Vec2D> nodePos);
    // Replaces nodes/edges/boundary data with a pre-built network, e.g. one read
    // from a Gmsh mesh. Gives every node of the supplied boundary chords tag 0.
    void initializeFromSpringMesh(const SpringMeshDescription& description);
    Index initializeBoundaryNodes(const BoundaryMeshBuildInformation& boundaryBuildInfo);
    void initializeInternalEdgesNodes(const InternalEdgesMeshBuildInformation internalEdgeBuildInfo);
    void addInternalNodes_Random(Index numNodesToAdd);
    void addInternalNodes_Uniform(int numNodesToTargetAdd);
    double calculateAreaOfMeshObject(Index numSamples);
    void createEdges(Index K, double distanceTheshold);

    // ---- boundary tagging (FVM identifyBoundaryEdgesAndCells analogue) --------
    // Returns a fresh, empty tag id. One per boundary polygon.
    Index createBoundaryTag();
    // Retags nodes inside `region`, always overriding (a node has one tag).
    // NO_BOUNDARY creates a new tag. Interior nodes are only promoted when
    // includeInteriorNodes is set.
    Index tagNodesByDomain(const Domain& region, Index boundaryId = NO_BOUNDARY,
                           bool includeInteriorNodes = false);
    // Adjacency, edge classification and normals. Idempotent.
    void rebuildBoundaryBookkeeping();
    void ensureBoundaryBookkeeping();
    void buildNodeEdgeAdjacency();
    void classifyBoundaryEdges();
    void computeBoundaryNodeNormals();

    // ---- boundary conditions --------------------------------------------------
    void setBoundaryCondition(Index boundaryId, const BoundaryCondition& condition);
    [[nodiscard]] const BoundaryCondition& boundaryCondition(Index boundaryId) const;
    [[nodiscard]] const std::vector<Index>& boundaryNodes(Index boundaryId) const;
    [[nodiscard]] bool isBoundaryNode(Index nodeIndex) const;
    void logBoundarySummary() const;

    // ---- flat vector (Eigen) --------------------------------------------------
    void createFlatVec();
    void cacheEdgesExtensionFromFlatVec();
    void syncVecWithFlatVecData();

private:
    struct EdgeKey
    {
        Index a;
        Index b;

        EdgeKey(Index i, Index j)
        {
            a = std::min(i, j);
            b = std::max(i, j);
        }

        bool operator==(const EdgeKey& other) const
        {
            return a == other.a && b == other.b;
        }
    };
    struct EdgeHash
    {
        size_t operator()(const EdgeKey& e) const
        {
            return std::hash<Index>()(e.a) ^
                   (std::hash<Index>()(e.b) << 1);
        }
    };

    std::vector<std::pair<Vec2D, Vec2D>> segmentDiscontinuities;
    Index numSegmentDiscontinuities = 0;
    bool areaCalculationDone = false;
    bool bookkeepingDirty = true;

    // Owned copies of the generation filters, with domainFilters pointing at them.
    std::vector<std::unique_ptr<Domain>> ownedDomainFilters;
    std::vector<const Domain*> domainFilters;
    std::vector<bool> domainFilterParity;

    void findNearestNeighbours(const KDTree& tree, const Vec2D& queryPoint, Index K, std::vector<size_t>& indices, std::vector<double>& distances);
    bool doesEdgeIntersectDiscontinuity(Index node1Idx, Index node2Idx);
    void resizeNodeAttributes();
    void rebuildBoundaryIdentifiers();
    void pruneBoundarySegments();
    // Rebuilds chords from spring adjacency for tags with no generation polygon.
    void inferBoundarySegments(Index boundaryId);
    void rebuildNormalsForBoundary(Index boundaryId, const Vec2D& normalOverride);
};
