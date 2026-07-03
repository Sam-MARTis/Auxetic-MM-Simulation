#pragma once

#include <vector>
#include <tuple>
#include <functional>
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
        
        struct NodeProperties{
            double x;
            double y;
        };
        struct Nodes{
            unsigned int numNodes;
            std::vector<NodeProperties> nodesPos; // N
            Constraints constraints;
        };
        void initializeBoundaryNodes()

    public:
        Nodes nodes;    
        Edges edges;
        double constraintParameter = 0.0f;
        void enforceConstraints(double t);
        void cacheEdgesExtension();
        Mesh();
};