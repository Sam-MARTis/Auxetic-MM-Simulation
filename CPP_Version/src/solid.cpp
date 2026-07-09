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
        const double dx = mesh. nodePosFlat[2*a] - mesh.nodePosFlat[2*b];
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


void Solid::solve()
{
    const int N   = mesh.nodes.numNodes;
    const int DOF = 2 * N;
    const int M   = A.rows();
    const int total = DOF + M;

    Vec lambda = Vec::Zero(M);    // initialise lagrange multipliers

    Eigen::SparseLU<SpMat> solver;

    for (int iter = 0; iter < MAX_SOLVER_INTERNAL_ITER; iter++)
    {
        // rebuild K and g from current positions
        assembleGradientAndHessian();

        // constraint violation
        Vec violation = A * mesh.nodePosFlat - c;

        // assemble augmented system
        SpMat K_aug(total, total);
        K_aug.reserve(K.nonZeros() + 2 * A.nonZeros());

        // Fill K_aug with K, A^T, and A
        for (int j = 0; j < K.outerSize(); j++)
            for (SpMat::InnerIterator it(K, j); it; ++it)
                K_aug.insert(it.row(), it.col()) = it.value();

        for (int j = 0; j < A.outerSize(); j++)
            for (SpMat::InnerIterator it(A, j); it; ++it)
            {
                K_aug.insert(it.col(), DOF + it.row()) = it.value();  // A^T
                K_aug.insert(DOF + it.row(), it.col()) = it.value();  // A
            }

        // RHS
        Vec rhs(total);
        rhs.head(DOF) = -g;
        rhs.tail(M)   = -violation;

        // factorize and solve
        solver.analyzePattern(K_aug);
        solver.factorize(K_aug);
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Factorization failed at iteration " << iter << "\n";
            break;
        }

        Vec sol = solver.solve(rhs);

        // update
        mesh.nodePosFlat += sol.head(DOF);
        lambda           += sol.tail(M);

        // convergence check
        double res_g = g.norm();
        double res_c = violation.norm();

        if (res_g < TOL && res_c < TOL)
        {
            std::cout << "Converged at iteration " << iter << "\n";
            break;
        }
    }
}

Solid::Solid(Mesh& solidMesh): mesh(solidMesh)
{
    std::vector<BoundaryEdgeInformation>& edgeInfo = mesh.boundaryEdgeInfo;
    std::vector<Vec2D> nodesPos = mesh.nodes.nodesPos;
    assert(!edgeInfo.empty());
    constraints.reserve(edgeInfo.size());
    for (unsigned int i=0; i<edgeInfo.size(); i++){
        assert(i == edgeInfo.at(i).id);
        const Vec2D& node1Pos = nodesPos.at(std::get<0>(edgeInfo.at(i).boundsIndices));
        const Vec2D& node2Pos = nodesPos.at(std::get<1>(edgeInfo.at(i).boundsIndices)-1);
        const Vec2D edgeVector = node2Pos - node1Pos;
        const Vec2D normalVec = Vec2D(edgeVector.y, -edgeVector.x).normalized();
        const double thisC = normalVec.x * node1Pos.x + normalVec.y * node1Pos.y;
        assert(std::abs(normalVec.x * node2Pos.x + normalVec.y * node2Pos.y - thisC) < TOL);
        BoundaryConstraint constraint;
        constraint.c = thisC;
        constraint.boundaryId = edgeInfo.at(i).id;
        constraint.nVec = normalVec;
        constraints.emplace_back(constraint);
    }
    mesh.createFlatVec();

    assembleConstraintMatrix();
    c = A * mesh.nodePosFlat;

}

void Solid::timeStep(const double dt)
{
    const float dC = dt*VELOCITY_CONSTRAINTS_EVOLUTION;
    assert(c.size() == mesh.numBoundaryPoints);
    for (unsigned int i = 0; i < mesh.numBoundaryPoints; i++){
        c[i] += dC;
    }

}

