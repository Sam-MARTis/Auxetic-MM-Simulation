#include "solid.hpp"

#include <iostream>
#include <omp.h>

void Solid::classifyDegreesOfFreedom()
{
    PROFILE_SCOPE("Solid::classifyDegreesOfFreedom");
    const Index numNodes = mesh.nodes.numNodes;
    const int numDofs = 2 * static_cast<int>(numNodes);

    freeDofToReduced.assign(numDofs, -1);
    reducedToFreeDof.clear();
    reducedToFreeDof.reserve(numDofs);
    numFixedNodes = 0;

    for (Index nodeIdx = 0; nodeIdx < numNodes; nodeIdx++)
    {
        const Index boundaryId = mesh.nodes.boundaryId.at(nodeIdx);
        const bool isFixed =
            boundaryId != NO_BOUNDARY &&
            mesh.boundaryConditions.at(boundaryId).type == BoundaryConditionType::FIXED;
        if (isFixed)
        {
            numFixedNodes++;
            continue;
        }
        for (int component = 0; component < 2; component++)
        {
            freeDofToReduced.at(2 * nodeIdx + component) = static_cast<int>(reducedToFreeDof.size());
            reducedToFreeDof.push_back(2 * nodeIdx + component);
        }
    }
}

void Solid::buildNormalConstraints()
{
    PROFILE_SCOPE("Solid::buildNormalConstraints");
    normalConstraints.clear();
    for (Index nodeIdx = 0; nodeIdx < mesh.nodes.numNodes; nodeIdx++)
    {
        const Index boundaryId = mesh.nodes.boundaryId.at(nodeIdx);
        if (boundaryId == NO_BOUNDARY) continue;
        const BoundaryCondition& condition = mesh.boundaryConditions.at(boundaryId);
        if (!condition.usesNormalConstraint()) continue;

        // One row per incident boundary plane: a straight edge contributes one, a
        // corner contributes one per edge, which pins it at the intersection of the
        // offset planes instead of letting it slide.
        for (Vec2D normal : mesh.nodes.constraintNormals.at(nodeIdx))
        {
            const double magnitude = normal.norm();
            if (magnitude <= EPSILON)
            {
                std::cerr << "Solid: node " << nodeIdx << " on boundary " << boundaryId
                          << " has no usable outward normal, skipping its plane constraint. "
                          << "Set BoundaryCondition::normalOverride for this boundary." << std::endl;
                continue;
            }
            normal /= magnitude;
            normalConstraints.push_back(
                NormalConstraint{boundaryId, nodeIdx, normal, condition.velocity});
        }
    }
}

void Solid::assembleConstraintMatrix()
{
    PROFILE_SCOPE("Solid::assembleConstraintMatrix");
    const int numConstraints = static_cast<int>(normalConstraints.size());
    const int numFreeDofs = static_cast<int>(reducedToFreeDof.size());

    std::vector<Eigen::Triplet<double>> constraintTriplets;
    constraintTriplets.reserve(static_cast<std::size_t>(numConstraints) * 2);
    for (int row = 0; row < numConstraints; row++)
    {
        const NormalConstraint& constraint = normalConstraints.at(row);
        const int xDof = freeDofToReduced.at(2 * constraint.nodeIndex);
        const int yDof = freeDofToReduced.at(2 * constraint.nodeIndex + 1);
        // Guard: an overlapping tag must not silently corrupt A.
        if (xDof < 0 || yDof < 0)
        {
            std::cerr << "Solid: node " << constraint.nodeIndex
                      << " is both fixed and plane-constrained; the plane constraint is ignored."
                      << std::endl;
            continue;
        }
        constraintTriplets.emplace_back(row, xDof, constraint.normal[0]);
        constraintTriplets.emplace_back(row, yDof, constraint.normal[1]);
    }
    A = SpMat(numConstraints, numFreeDofs);
    A.setFromTriplets(constraintTriplets.begin(), constraintTriplets.end());
}

void Solid::buildConstraintRHS()
{
    PROFILE_SCOPE("Solid::buildConstraintRHS");
    const int numConstraints = static_cast<int>(normalConstraints.size());
    c = VecXd(numConstraints);
    for (int row = 0; row < numConstraints; row++)
    {
        const NormalConstraint& constraint = normalConstraints.at(row);
        // Built from the live flat vector, so refreshBoundaryConditions() works.
        const Vec2D position{mesh.nodePosFlat[2 * constraint.nodeIndex],
                             mesh.nodePosFlat[2 * constraint.nodeIndex + 1]};
        c[row] = Vec2D::dot(constraint.normal, position);
    }
}

VecXd Solid::gatherReducedPositions() const
{
    VecXd reduced(reducedToFreeDof.size());
    for (std::size_t i = 0; i < reducedToFreeDof.size(); i++)
    {
        reduced[static_cast<int>(i)] = mesh.nodePosFlat[reducedToFreeDof[i]];
    }
    return reduced;
}

void Solid::assembleGradientAndHessian()
{
    PROFILE_SCOPE("Solid::assembleGradientAndHessian");
    mesh.cacheEdgesExtensionFromFlatVec();
    const int numFreeDofs = static_cast<int>(reducedToFreeDof.size());

    g = VecXd::Zero(numFreeDofs);

    const int nThreads = omp_get_max_threads();
    std::vector<VecXd> g_local(nThreads, VecXd::Zero(numFreeDofs));
    std::vector<std::vector<Eigen::Triplet<double>>> triplets_local(nThreads);
    for (int t = 0; t < nThreads; t++)
        triplets_local[t].reserve((mesh.edges.numEdges / nThreads + 1) * 16);

    const std::vector<EdgesConnection>& connections = mesh.edges.edgesConnection;
    const std::vector<EdgeProperties>& props = mesh.edges.edgesProperty;
    const std::vector<int>& dofMap = freeDofToReduced;

    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < mesh.edges.numEdges; i++)
    {
        const int tid = omp_get_thread_num();

        const int a = static_cast<int>(connections[i].node1);
        const int b = static_cast<int>(connections[i].node2);
        const double dx = mesh.nodePosFlat[2*a]   - mesh.nodePosFlat[2*b];
        const double dy = mesh.nodePosFlat[2*a+1] - mesh.nodePosFlat[2*b+1];
        const double d  = props[i].cached_extension;
        const double l0 = props[i].natural_length;
        const double k  = props[i].stiffness;

        if (d < EPSILON) continue;

        const double inv   = 1.0 / d;
        const double alpha = k * (1.0 - l0 * inv);
        const double beta  = k * l0 * inv;
        const double fx    = alpha * dx;
        const double fy    = alpha * dy;

        // Only free DOFs get gradient rows; a fixed neighbour still pushes them.
        const int g1x = 2*a,   g1y = 2*a+1;
        const int g2x = 2*b,   g2y = 2*b+1;
        const int r1x = dofMap[g1x], r1y = dofMap[g1y];
        const int r2x = dofMap[g2x], r2y = dofMap[g2y];
        if (r1x >= 0) g_local[tid][r1x] += fx;
        if (r1y >= 0) g_local[tid][r1y] += fy;
        if (r2x >= 0) g_local[tid][r2x] -= fx;
        if (r2y >= 0) g_local[tid][r2y] -= fy;

        const double ex  = dx * inv;
        const double ey  = dy * inv;
        const double c00 = alpha + beta * ex * ex;
        const double c01 = beta  * ex * ey;
        const double c11 = alpha + beta * ey * ey;

        auto& tl = triplets_local[tid];
        // Takes *global* DOF indices; do not pre-reduce (DeveloperReference 4.4).
        const auto addEntry = [&tl, &dofMap](const int row, const int col, const double value)
        {
            const int reducedRow = dofMap[row];
            const int reducedCol = dofMap[col];
            if (reducedRow < 0 || reducedCol < 0) return;
            tl.emplace_back(reducedRow, reducedCol, value);
        };

        addEntry(g1x, g1x,  c00);
        addEntry(g1x, g1y,  c01);
        addEntry(g1y, g1x,  c01);
        addEntry(g1y, g1y,  c11);
        addEntry(g2x, g2x,  c00);
        addEntry(g2x, g2y,  c01);
        addEntry(g2y, g2x,  c01);
        addEntry(g2y, g2y,  c11);
        addEntry(g1x, g2x, -c00);
        addEntry(g1x, g2y, -c01);
        addEntry(g1y, g2x, -c01);
        addEntry(g1y, g2y, -c11);
        addEntry(g2x, g1x, -c00);
        addEntry(g2x, g1y, -c01);
        addEntry(g2y, g1x, -c01);
        addEntry(g2y, g1y, -c11);
    }

    for (int t = 0; t < nThreads; t++)
        g += g_local[t];

    triplets.clear();
    for (int t = 0; t < nThreads; t++)
        triplets.insert(triplets.end(), triplets_local[t].begin(), triplets_local[t].end());

    K = SpMat(numFreeDofs, numFreeDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());
}

void Solid::applyBoundaryVelocity(const double dt)
{
    for (std::size_t row = 0; row < normalConstraints.size(); row++)
    {
        c[static_cast<int>(row)] += normalConstraints.at(row).velocity * dt;
    }
}

void Solid::solve()
{
    PROFILE_SCOPE("Solid::solve");
    const int numFreeDofs = static_cast<int>(reducedToFreeDof.size());
    const int numConstraints = static_cast<int>(normalConstraints.size());
    const int total = numFreeDofs + numConstraints;
    if (total == 0) return;

    VecXd lambda = VecXd::Zero(numConstraints);

    for (int iter = 0; iter < MAX_SOLVER_INTERNAL_ITER; iter++)
    {
        {
            PROFILE_SCOPE("Solid::solve / assembleGradientAndHessian");
            assembleGradientAndHessian();
        }

        VecXd violation = VecXd::Zero(numConstraints);
        {
            PROFILE_SCOPE("Solid::solve / constraintViolation");
            if (numConstraints > 0) violation = A * gatherReducedPositions() - c;
        }

        {
            PROFILE_SCOPE("Solid::solve / assembleAugmentedSystem");
            std::vector<Eigen::Triplet<double>> augmentedTriplets;
            augmentedTriplets.reserve(K.nonZeros() + 2 * A.nonZeros());

            for (int j = 0; j < K.outerSize(); j++)
                for (SpMat::InnerIterator it(K, j); it; ++it)
                    augmentedTriplets.emplace_back(it.row(), it.col(), it.value());

            for (int j = 0; j < A.outerSize(); j++)
                for (SpMat::InnerIterator it(A, j); it; ++it)
                {
                    augmentedTriplets.emplace_back(it.col(), numFreeDofs + it.row(), it.value());
                    augmentedTriplets.emplace_back(numFreeDofs + it.row(), it.col(), it.value());
                }

            K_aug = SpMat(total, total);
            K_aug.setFromTriplets(augmentedTriplets.begin(), augmentedTriplets.end());
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
        rhs.head(numFreeDofs) = -g;
        if (numConstraints > 0) rhs.tail(numConstraints) = -violation;

        VecXd sol;
        {
            PROFILE_SCOPE("Solid::solve / solveLinearSystem");
            sol = solver.solve(rhs);
        }

        {
            PROFILE_SCOPE("Solid::solve / updateState");
            // Fixed DOFs are absent from the map, so they stay pinned.
            for (int i = 0; i < numFreeDofs; i++)
                mesh.nodePosFlat[reducedToFreeDof[i]] += sol[i];
            if (numConstraints > 0) lambda += sol.tail(numConstraints);
        }

        const double res_g = g.norm();
        const double res_c = numConstraints > 0 ? violation.norm() : 0.0;

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
    if (mesh.boundaryConditions.size() != mesh.numBoundaries)
    {
        mesh.boundaryConditions.resize(mesh.numBoundaries, BoundaryCondition::free());
    }
    mesh.ensureBoundaryBookkeeping();
    mesh.createFlatVec();
    refreshBoundaryConditions();

    triplets.reserve(static_cast<std::size_t>(mesh.edges.numEdges) * 16);

    std::cout << "Solid: " << mesh.nodes.numNodes << " nodes ("
              << numFixedNodes << " fixed, " << normalConstraints.size()
              << " plane-constrained, "
              << reducedToFreeDof.size() << " free DOFs)" << std::endl;
}

void Solid::timeStep(const double dt)
{
    PROFILE_SCOPE("Solid::timeStep");
    applyBoundaryVelocity(dt);
    solve();
}

void Solid::refreshBoundaryConditions()
{
    PROFILE_SCOPE("Solid::refreshBoundaryConditions");
    mesh.ensureBoundaryBookkeeping();
    classifyDegreesOfFreedom();
    buildNormalConstraints();
    assembleConstraintMatrix();
    buildConstraintRHS();
    // Reduced size and sparsity changed.
    patternAnalyzed = false;
}

double Solid::computeElasticEnergyFromExtensionCache()
{
    const Index numEdges = mesh.edges.numEdges;
    const std::vector<EdgeProperties>& edgeProps = mesh.edges.edgesProperty;
    double totalEnergy = 0.0;

    for (Index i = 0; i < numEdges; i++)
    {
        const double d  = edgeProps[i].cached_extension;
        const double l0 = edgeProps[i].natural_length;
        const double k  = edgeProps[i].stiffness;

        if (d < EPSILON) continue;

        totalEnergy += 0.5 * k * (d - l0) * (d - l0);
    }
    return totalEnergy;
}

void Solid::logCurrentElasticEnergy()
{
    const double energy = computeElasticEnergyFromExtensionCache();
    EnergyHistory.push_back(energy);
}

void Solid::logBoundaryConditionSummary() const
{
    std::cout << "Solver boundary treatment:" << std::endl;
    for (std::size_t row = 0; row < normalConstraints.size(); row++)
    {
        const NormalConstraint& constraint = normalConstraints.at(row);
        if (row < 4 || row + 4 >= normalConstraints.size())
        {
            std::cout << "  constraint " << row
                      << " | boundary " << constraint.boundaryId
                      << " | node " << constraint.nodeIndex
                      << " | n (" << constraint.normal[0] << ", " << constraint.normal[1] << ")"
                      << " | velocity " << constraint.velocity
                      << " | c0 " << c[static_cast<int>(row)] << std::endl;
        }
        else if (row == 4)
        {
            std::cout << "  ... (" << normalConstraints.size() - 8 << " more)" << std::endl;
        }
    }
    std::cout << "  fixed nodes: " << numFixedNodes
              << ", plane-constrained nodes: " << normalConstraints.size()
              << ", free DOFs: " << reducedToFreeDof.size() << std::endl;
}

double Solid::constraintViolation() const
{
    if (c.size() == 0 || A.rows() == 0) return 0.0;
    const VecXd residual = A * gatherReducedPositions() - c;
    return residual.cwiseAbs().maxCoeff();
}

double Solid::gradientNorm() const
{
    return g.size() == 0 ? 0.0 : g.norm();
}
