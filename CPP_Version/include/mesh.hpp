#pragma once


#include <tuple>
#include <functional>

#include "utils.hpp"
#include "loader.hpp"

#include <iostream>
#include <nanoflann.hpp>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/Dense>


struct PointCloud
{
    std::vector<Vec2D>& pointsPos;
    size_t kdtree_get_point_count() const { return pointsPos.size(); };
    double kdtree_get_pt(const size_t idx, int dim) const { return dim==0?pointsPos[idx].x:pointsPos[idx].y; };
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false;}
};
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>;

struct BoundaryEdgeInformation
{
    unsigned int id;
    std::pair<unsigned int, unsigned int> boundsIndices;
};

struct EdgesConnection{
    unsigned int node1;
    unsigned int node2;
};

struct EdgeProperties{
    double natural_length;
    double stiffness;
    double cached_extension;
};
class Mesh{
    private:
        struct ConstraintParameter{
            double x0;
            double m;
        };

        struct Edges{
            unsigned int numEdges = 0;
            std::vector<EdgeProperties> edgesProperty; // M
            std::vector<EdgesConnection> edgesConnection; // M
        };

        // struct Constraint{
        //     std::function<double(double x0t, double t, double m)> func;
        //     ConstraintParameter param;
        // };
        // struct Constraints{
        //     unsigned int numConstraints = 0;
        //     unsigned int numFree = 0;
        //     std::vector<unsigned int> fixedIndices; // numConstraints
        //     std::vector<unsigned int> freeIndices; // numNodes - numConstraints
        //     // std::vector<Constraint> constraints; // numConstraints
        // };


        struct Nodes{
            unsigned int numNodes = 0;
            std::vector<Vec2D> nodesPos; // N
            // Constraints constraints;
        };

        struct EdgeKey
        {
            unsigned int a;
            unsigned int b;

            EdgeKey(unsigned int i, unsigned int j)
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
                return std::hash<unsigned int>()(e.a) ^
                       (std::hash<unsigned int>()(e.b) << 1);
            }
        };
        bool areaCalculationDone = false;



        // nanoflann::KNNResultSet<double> result;
        // KDTree tree;
        std::vector<std::pair<Vec2D, Vec2D>> segmentDiscontinuities;
        unsigned int numSegmentDiscontinuities = 0;
        void findNearestNeighbours(const KDTree& tree, const Vec2D& queryPoint, unsigned int K, std::vector<size_t>& indices, std::vector<double>& distances);
        bool doesEdgeIntersectDiscontinuity(const unsigned int node1Idx, unsigned int node2Idx);
    public:
        Nodes nodes;
        Vec nodePosFlat;
        Edges edges;
        double areaOfDxf = 4.0;
        unsigned int highestBoundaryID = 0;
        std::vector<BoundaryEdgeInformation> boundaryEdgeInfo;
        unsigned int numBoundaries = 0;
        unsigned int numBoundaryPoints = 0;
        std::vector<Boundary> boundaries;
        std::vector<bool> boundaryIncludeExclude;
        void insertNodesAtPositions(const std::vector<Vec2D> nodePos);
        void initializeBoundaryNodes(const BoundaryMeshBuildInformation& boundaryBuildInfo);
        void initializeInternalEdgesNodes(const InternalEdgesMeshBuildInformation internalEdgeBuildInfo);
        void addInternalNodes_Random(const unsigned int numNodesToAdd);
        void addInternalNodes_Uniform(const int numNodesToTargetAdd);
        double calculateAreaOfMeshObject(const unsigned int numSamples);
        double constraintParameter = 0.0f;
        void createEdges(const unsigned int K, const double distanceTheshold);
        // void enforceConstraints(double t);
        void createFlatVec();
        void cacheEdgesExtensionFromFlatVec();
        void syncVecWithFlatVecData();


        // Mesh();
};



