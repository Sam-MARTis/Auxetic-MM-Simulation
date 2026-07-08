#pragma once

#include <vector>
#include <tuple>
#include <functional>

#include "utils.hpp"
#include "loader.hpp"

struct BoundaryEdgeInformation
{
    unsigned int id;
    std::pair<unsigned int, unsigned int> boundsIndices;
};

class Mesh{
    private:
        struct ConstraintParameter{
            double x0;
            double m;
        };
        struct EdgeProperties{
            double natural_length;
            double stiffness;
            double cached_extension;
        };
        struct EdgesConnection{
            unsigned int node1;
            unsigned int node2;
        };
        struct Edges{
            unsigned int numEdges;
            std::vector<EdgeProperties> edgesProperty; // M
            std::vector<EdgesConnection> edgesConnection; // M
        };

        struct Constraint{
            std::function<double(double x0t, double t, double m)> func;
            ConstraintParameter param;
        };
        struct Constraints{
            unsigned int numConstraints;
            unsigned int numFree;
            std::vector<unsigned int> fixedIndices; // numConstraints
            std::vector<unsigned int> freeIndices; // numNodes - numConstraints
            std::vector<Constraint> constraints; // numConstraints
        };
        

        struct Nodes{
            unsigned int numNodes;
            std::vector<Vec2D> nodesPos; // N
            Constraints constraints;
        };
        bool areaCalculationDone = false;

    public:
        Nodes nodes;
        Edges edges;
        double areaOfDxf = 4.0;
        unsigned int highestBoundaryID = 0;
        std::vector<BoundaryEdgeInformation> boundaryEdgeInfo;
        unsigned int numBoundaries = 0;
        std::vector<Boundary> boundaries;
        std::vector<bool> boundaryIncludeExclude;
        void insertNodesAtPositions(const std::vector<Vec2D> nodePos);
        void initializeBoundaryNodes(const BoundaryMeshBuildInformation& boundaryBuildInfo);
        void initializeInternalEdgesNodes(const InternalEdgesMeshBuildInformation internalEdgeBuildInfo);
        void addInternalNodes_Random(const unsigned int numNodesToAdd);
        void addInternalNodes_Uniform(const int numNodesToTargetAdd);
        double calculateAreaOfMeshObject(const unsigned int numSamples);
        double constraintParameter = 0.0f;
        void createEdges(const double distanceTheshold);
        void enforceConstraints(double t);
        void cacheEdgesExtension();

        // Mesh();
};

