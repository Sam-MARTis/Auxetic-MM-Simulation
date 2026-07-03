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
            std::vector<EdgeProperties> edges; // M
            std::vector<EdgesConnection> edges_connection; // M
        };

    public:
        unsigned int numNodes; // N
        std::vector<double> positions; // 2*N
        std::vector<unsigned int> freeIndices;
        std::vector<unsigned int> fixedIndices;
        // vector of functions with x0, t and m
        std::vector<std::function<double(double x0t, double t, double m)>> constraintFunctions; // fixedIndices num
        std::vector<ConstraintParameter> constraintParameters; 

        unsigned int numFree;
        unsigned int numFixed;

        // unsigned int numEdges; // M
        // std::vector<double> natural_lengths; // M
        // std::vector<double> cached_extensions; // M
        // std::vector<double> stiffnesses; // M

        // std::vector<std::tuple<unsigned int, unsigned int>> edges_connection; // M        
        Edges edges;

        void enforceConstraints(std::vector<double> &positions, double t);
        Mesh();

};