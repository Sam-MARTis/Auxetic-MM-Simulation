#include "solid.hpp"


void Solid::assembleConstraintMatrix()
{
    PROFILE_SCOPE("Solid::assembleConstraintMatrix");
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
            triplets.emplace_back(row, 2*i, con.nVec[0]);
            triplets.emplace_back(row, 2*i+1, con.nVec[1]);
            assert(row < (int)M);
            assert(2*i+1 < 2*mesh.nodes.numNodes);
            row++;
        }
    }
    A = SpMat(M, 2*mesh.nodes.numNodes);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void Solid::buildConstraintRHS()
{
    PROFILE_SCOPE("Solid::buildConstraintRHS");
    const unsigned M = mesh.numBoundaryPoints;
    c = VecXd(M);
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
    PROFILE_SCOPE("Solid::assembleGradientAndHessian");
    mesh.cacheEdgesExtensionFromFlatVec();
    const int N   = mesh.nodes.numNodes;
    const int DOF = 2 * N;

    g = VecXd::Zero(DOF);

    const int nThreads = omp_get_max_threads();
    std::vector<VecXd> g_local(nThreads, VecXd::Zero(DOF));
    std::vector<std::vector<Eigen::Triplet<double>>> triplets_local(nThreads);
    for (int t = 0; t < nThreads; t++)
        triplets_local[t].reserve((mesh.edges.numEdges / nThreads + 1) * 16);

    std::vector<EdgesConnection>& connections = mesh.edges.edgesConnection;
    std::vector<EdgeProperties>&  props       = mesh.edges.edgesProperty;

    #pragma omp parallel for schedule(static)
    for (unsigned int i = 0; i < mesh.edges.numEdges; i++)
    {
        const int tid = omp_get_thread_num();

        const int a  = connections[i].node1;
        const int b  = connections[i].node2;
        const double dx = mesh.nodePosFlat[2*a]   - mesh.nodePosFlat[2*b];
        const double dy = mesh.nodePosFlat[2*a+1] - mesh.nodePosFlat[2*b+1];
        const double d  = props[i].cached_extension;   // keep cache
        const double l0 = props[i].natural_length;
        const double k  = props[i].stiffness;

        if (d < EPSILON) continue;

        const double inv   = 1.0 / d;
        const double alpha = k * (1.0 - l0 * inv);
        const double beta  = k * l0 * inv;
        const double fx    = alpha * dx;
        const double fy    = alpha * dy;

        g_local[tid][2*a]   += fx;
        g_local[tid][2*a+1] += fy;
        g_local[tid][2*b]   -= fx;
        g_local[tid][2*b+1] -= fy;

        const double ex  = dx * inv;
        const double ey  = dy * inv;
        const double c00 = alpha + beta * ex * ex;
        const double c01 = beta  * ex * ey;
        const double c11 = alpha + beta * ey * ey;

        const int r1x = 2*a, r1y = 2*a+1;
        const int r2x = 2*b, r2y = 2*b+1;

        auto& tl = triplets_local[tid];
        tl.emplace_back(r1x, r1x,  c00);
        tl.emplace_back(r1x, r1y,  c01);
        tl.emplace_back(r1y, r1x,  c01);
        tl.emplace_back(r1y, r1y,  c11);
        tl.emplace_back(r2x, r2x,  c00);
        tl.emplace_back(r2x, r2y,  c01);
        tl.emplace_back(r2y, r2x,  c01);
        tl.emplace_back(r2y, r2y,  c11);
        tl.emplace_back(r1x, r2x, -c00);
        tl.emplace_back(r1x, r2y, -c01);
        tl.emplace_back(r1y, r2x, -c01);
        tl.emplace_back(r1y, r2y, -c11);
        tl.emplace_back(r2x, r1x, -c00);
        tl.emplace_back(r2x, r1y, -c01);
        tl.emplace_back(r2y, r1x, -c01);
        tl.emplace_back(r2y, r1y, -c11);
    }

    // merge g
    for (int t = 0; t < nThreads; t++)
        g += g_local[t];

    // merge and build K
    triplets.clear();
    for (int t = 0; t < nThreads; t++)
        triplets.insert(triplets.end(),
                        triplets_local[t].begin(),
                        triplets_local[t].end());

    K = SpMat(DOF, DOF);
    K.setFromTriplets(triplets.begin(), triplets.end());
}

void Solid::solve()
{
    PROFILE_SCOPE("Solid::solve");
    const int N     = mesh.nodes.numNodes;
    const int DOF   = 2 * N;
    const int M     = A.rows();
    const int total = DOF + M;

    VecXd lambda = VecXd::Zero(M);

    for (int iter = 0; iter < MAX_SOLVER_INTERNAL_ITER; iter++)
    {
        {
            PROFILE_SCOPE("Solid::solve / assembleGradientAndHessian");
            assembleGradientAndHessian();
        }

        VecXd violation;
        {
            PROFILE_SCOPE("Solid::solve / constraintViolation");
            violation = A * mesh.nodePosFlat - c;
        }

        {
            PROFILE_SCOPE("Solid::solve / assembleAugmentedSystem");
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(K.nonZeros() + 2 * A.nonZeros());

            for (int j = 0; j < K.outerSize(); j++)
                for (SpMat::InnerIterator it(K, j); it; ++it)
                    triplets.emplace_back(it.row(), it.col(), it.value());

            for (int j = 0; j < A.outerSize(); j++)
                for (SpMat::InnerIterator it(A, j); it; ++it)
                {
                    triplets.emplace_back(it.col(), DOF + it.row(), it.value());
                    triplets.emplace_back(DOF + it.row(), it.col(), it.value());
                }

            K_aug = SpMat(total, total);
            K_aug.setFromTriplets(triplets.begin(), triplets.end());
        }

        {
            PROFILE_SCOPE("Solid::solve / analyzePattern");
            if (!patternAnalyzed)
            {
                solver.analyzePattern(K_aug);
                patternAnalyzed = true;
            }
        }

        {
            PROFILE_SCOPE("Solid::solve / factorize");
            solver.factorize(K_aug);
        }

        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Factorization failed at iteration " << iter << "\n";
            break;
        }

        VecXd rhs(total);
        rhs.head(DOF) = -g;
        rhs.tail(M)   = -violation;

        VecXd sol;
        {
            PROFILE_SCOPE("Solid::solve / solveLinearSystem");
            sol = solver.solve(rhs);
        }

        {
            PROFILE_SCOPE("Solid::solve / updateState");
            mesh.nodePosFlat += sol.head(DOF);
            lambda           += sol.tail(M);
        }

        double res_g, res_c;
        {
            PROFILE_SCOPE("Solid::solve / convergenceCheck");
            res_g = g.norm();
            res_c = violation.norm();
        }

        if (res_g < TOL && res_c < TOL)
        {
            std::cout << "Converged at iteration " << iter << "\n";
            break;
        }
    }
}

Solid::Solid(Mesh& solidMesh): mesh(solidMesh)
{
    PROFILE_SCOPE("Solid::Solid");
    std::cout << "Solid1"<<std::endl;
    std::vector<BoundaryEdgeInformation>& edgeInfo = mesh.boundaryEdgeInfo;
    std::vector<Vec2D> nodesPos = mesh.nodes.nodesPos;
    assert(!edgeInfo.empty());
    constraints.reserve(edgeInfo.size());
    std::cout << "Solid2"<<std::endl;
    for (unsigned int i=0; i<edgeInfo.size(); i++){
        assert(i == edgeInfo.at(i).id);
        const Vec2D& node1Pos = nodesPos.at(std::get<0>(edgeInfo.at(i).boundsIndices));
        const Vec2D& node2Pos = nodesPos.at(std::get<1>(edgeInfo.at(i).boundsIndices)-1);
        const Vec2D edgeVector = node2Pos - node1Pos;
        const Vec2D normalVec = Vec2D{-edgeVector[1], edgeVector[0]}.normalized();
        const double thisC = normalVec[0] * node1Pos[0] + normalVec[1] * node1Pos[1];
        assert(std::abs(normalVec[0] * node2Pos[0] + normalVec[1] * node2Pos[1] - thisC) < TOL);
        BoundaryConstraint constraint;
        constraint.c = thisC;
        constraint.boundaryId = edgeInfo.at(i).id;
        constraint.nVec = normalVec;
        constraints.emplace_back(constraint);
    }
    std::cout << "Solid3"<<std::endl;
    mesh.createFlatVec();
    std::cout << "nodePosFlat size: " << mesh.nodePosFlat.size() << "\n";
    std::cout << "expected: " << 2 * mesh.nodes.numNodes << "\n";
    assembleConstraintMatrix();
    c = A * mesh.nodePosFlat;
    std::cout<<"c[0] initial : "<<c[0]<<std::endl;
    triplets.reserve(mesh.edges.numEdges * 16);

}

void Solid::timeStep(const double dt)
{
    PROFILE_SCOPE("Solid::timeStep");
    const float dC = dt*VELOCITY_CONSTRAINTS_EVOLUTION;
    assert(c.size() == mesh.numBoundaryPoints);
    for (unsigned int i = 0; i < mesh.numBoundaryPoints; i++){
        c[i] += dC;
    }
    solve();
    // std::cout<<"c[0] : "<<c[0]<<std::endl;
    

}

double Solid::computeElasticEnergyFromExtensionCache()
{
    const unsigned int numEdges = mesh.edges.numEdges;
    const std::vector<EdgeProperties>& edgeProps = mesh.edges.edgesProperty;
    double totalEnergy = 0.0;
    
    for (unsigned int i = 0; i < numEdges; i++)
    {
        const double d  = edgeProps[i].cached_extension;
        const double l0 = edgeProps[i].natural_length;
        const double k  = edgeProps[i].stiffness;

        if (d < EPSILON) continue;

        const double energy = 0.5 * k * (d - l0) * (d - l0);
        totalEnergy += energy;
    }
    return totalEnergy;
}

void Solid::logCurrentElasticEnergy()
{
    const double energy = computeElasticEnergyFromExtensionCache();
    EnergyHistory.push_back(energy);
}