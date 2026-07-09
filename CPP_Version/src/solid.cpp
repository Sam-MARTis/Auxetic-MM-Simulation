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
