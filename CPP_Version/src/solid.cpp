#include "solid.hpp"


void Solid::assembleConstraintMatrix()
{
    const unsigned int M = mesh.numBoundaryPoints;
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(M * 2);
    int row = 0;
    for (const BoundaryConstraint& con : constraints)
    {
        const BoundaryEdgeInformation& boundary = mesh.boundaryEdgeInfo.at(con.boundaryId);
        const unsigned int startingIdx = std::get<0>(boundary.boundsIndices);
        const unsigned int endingIdx = std::get<1>(boundary.boundsIndices);

        for (unsigned int i=startingIdx; i<endingIdx; i++)
        {
            triplets.emplace_back(row, 2*i, con.nVec.x);
            triplets.emplace_back(row, 2*i+1, con.nVec.y);
            row++;
        }
    }
    A = SpMat(M, 2*mesh.nodes.numNodes);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void Solid::buildConstraintRHS()
{
    const unsigned M = mesh.numBoundaryPoints;
    c = Vec(M);
    int row = 0;
    for (const BoundaryConstraint& con : constraints)
    {
        const BoundaryEdgeInformation& boundary = mesh.boundaryEdgeInfo.at(con.boundaryId);
        const unsigned int startingIdx = std::get<0>(boundary.boundsIndices);
        const unsigned int endingIdx = std::get<1>(boundary.boundsIndices);
        unsigned int count = endingIdx - startingIdx;
        c.segment(row, count).setConstant(con.c);
        row += count;
    }
}

void Solid::assembleGradientAndHessian()
{
    mesh.cacheEdgesExtensionFromFlatVec();
    const int N = mesh.nodes.numNodes;
    const int DOF = 2*N;
    g = Vec::Zero(DOF);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.edges.numEdges * 16);
    std::vector<EdgesConnection>& connections = mesh.edges.edgesConnection;
    std::vector<EdgeProperties>& props = mesh.edges.edgesProperty;
    assert(connections.size() == mesh.edges.numEdges);
   
    for (unsigned int i=0; i< mesh.edges.numEdges; i++){
        const EdgesConnection& thisEdge = connections.at(i);
        const EdgeProperties& thisProp = props.at(i);
        const int a = thisEdge.node1;
        const int b = thisEdge.node2;
        const double dx = mesh.nodePosFlat[2*a] - mesh.nodePosFlat[2*b];
        const double dy = mesh.nodePosFlat[2*a+1] - mesh.nodePosFlat[2*b+1];
        const double d = thisProp.cached_extension;
        const double l0 = thisProp.natural_length;
        const double k = thisProp.stiffness;
        if (d < EPSILON) continue;
        const double inv = 1.0/d;
        const double alpha = k * (1.0 - l0 * inv);
        const double beta = k * l0 * inv;
        const double fx = alpha * dx;
        const double fy = alpha * dy;
        g[2*a] += fx;
        g[2*a+1] += fy;
        g[2*b] -= fx;
        g[2*b+1] -= fy;
        const double ex = dx * inv;
        const double ey = dy * inv;
        const double c00 = alpha + beta * ex * ex;
        const double c01 = beta * ex * ey;
        const double c11 = alpha + beta * ey * ey;
        const int r1x = 2*a, r1y = 2*a+1;
        const int r2x = 2*b, r2y = 2*b+1;
        triplets.emplace_back(r1x, r1x, c00);
        triplets.emplace_back(r1x, r1y, c01);
        triplets.emplace_back(r1y, r1x, c01);
        triplets.emplace_back(r1y, r1y, c11);

        triplets.emplace_back(r2x, r2x, c00);
        triplets.emplace_back(r2x, r2y, c01);
        triplets.emplace_back(r2y, r2x, c01);
        triplets.emplace_back(r2y, r2y, c11);

        triplets.emplace_back(r1x, r2x, -c00);
        triplets.emplace_back(r1x, r2y, -c01);
        triplets.emplace_back(r1y, r2x, -c01);
        triplets.emplace_back(r1y, r2y, -c11);

        triplets.emplace_back(r2x, r1x, -c00);
        triplets.emplace_back(r2x, r1y, -c01);
        triplets.emplace_back(r2y, r1x, -c01);
        triplets.emplace_back(r2y, r1y, -c11);



    }
    K = SpMat(DOF, DOF);
    K.setFromTriplets(triplets.begin(), triplets.end());
}
