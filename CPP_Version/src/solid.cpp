#include "solid.hpp"

#include <cmath>
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
    mesh.cacheTorsionStateFromFlatVec();

    if (mesh.numNearCollinearTorsions > 0)
    {
        if (!nearCollinearWarned)
        {
            std::cerr << "Solid: " << mesh.numNearCollinearTorsions
                      << " torsional element(s) are near-collinear; the cosine "
                         "energy is stationary there. Check the mesh/deformation."
                      << std::endl;
            nearCollinearWarned = true;
        }
    }
    else
    {
        nearCollinearWarned = false;
    }

    const int numFreeDofs = static_cast<int>(reducedToFreeDof.size());
    g = VecXd::Zero(numFreeDofs);

    const int nThreads = omp_get_max_threads();
    std::vector<VecXd> g_local(nThreads, VecXd::Zero(numFreeDofs));
    std::vector<std::vector<Eigen::Triplet<double>>> triplets_local(nThreads);
    const std::size_t tripletsPerThread =
        (static_cast<std::size_t>(mesh.edges.numEdges) * 16
         + static_cast<std::size_t>(mesh.torsions.numTorsions) * 36) / nThreads + 1;
    for (int t = 0; t < nThreads; t++)
        triplets_local[t].reserve(tripletsPerThread);

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

    // ---- bond-angle (torsional) elements -------------------------------------
    const std::vector<TorsionConnection>& torsionConnections = mesh.torsions.torsionConnection;
    const std::vector<TorsionProperties>& torsionProps = mesh.torsions.torsionProperty;

    #pragma omp parallel for schedule(static)
    for (Index t = 0; t < mesh.torsions.numTorsions; t++)
    {
        const int tid = omp_get_thread_num();
        const TorsionConnection& connection = torsionConnections[t];
        const TorsionProperties& property = torsionProps[t];

        const double a = property.cached_a;
        const double b = property.cached_b;
        if (a < EPSILON || b < EPSILON) continue;

        const int i = static_cast<int>(connection.nodeI);
        const int j = static_cast<int>(connection.nodeApex);
        const int k = static_cast<int>(connection.nodeK);

        const double ux = mesh.nodePosFlat[2*i]   - mesh.nodePosFlat[2*j];
        const double uy = mesh.nodePosFlat[2*i+1] - mesh.nodePosFlat[2*j+1];
        const double vx = mesh.nodePosFlat[2*k]   - mesh.nodePosFlat[2*j];
        const double vy = mesh.nodePosFlat[2*k+1] - mesh.nodePosFlat[2*j+1];

        const double ex = ux / a;
        const double ey = uy / a;
        const double fx = vx / b;
        const double fy = vy / b;

        const double c = property.cached_cosine;
        const double s = property.cached_sine;
        const double kappa = property.stiffness;
        const double c0 = property.natural_cosine;
        const double U1 = kappa * (c - c0);
        const double U2 = kappa;

        // J e = (-ey, ex), J f = (-fy, fx)
        const double Jex = -ey, Jey = ex;
        const double Jfx = -fy, Jfy = fx;
        const double gux = (s / a) * Jex, guy = (s / a) * Jey;
        const double gvx = -(s / b) * Jfx, gvy = -(s / b) * Jfy;

        const int dix = 2*i, diy = 2*i+1;
        const int djx = 2*j, djy = 2*j+1;
        const int dkx = 2*k, dky = 2*k+1;
        const int rix = dofMap[dix], riy = dofMap[diy];
        const int rjx = dofMap[djx], rjy = dofMap[djy];
        const int rkx = dofMap[dkx], rky = dofMap[dky];
        if (rix >= 0) g_local[tid][rix] += U1 * gux;
        if (riy >= 0) g_local[tid][riy] += U1 * guy;
        if (rjx >= 0) g_local[tid][rjx] -= U1 * (gux + gvx);
        if (rjy >= 0) g_local[tid][rjy] -= U1 * (guy + gvy);
        if (rkx >= 0) g_local[tid][rkx] += U1 * gvx;
        if (rky >= 0) g_local[tid][rky] += U1 * gvy;

        // Cosine Hessian blocks, physicsDerivation.md section 5.3.
        const double invA2 = 1.0 / (a * a);
        const double invB2 = 1.0 / (b * b);
        const double P00 = -(2.0*ex*fx + c*(1.0 - 3.0*ex*ex)) * invA2;
        const double P01 = -(ex*fy + fx*ey - 3.0*c*ex*ey) * invA2;
        const double P11 = -(2.0*ey*fy + c*(1.0 - 3.0*ey*ey)) * invA2;

        const double Q00 = -(2.0*fx*ex + c*(1.0 - 3.0*fx*fx)) * invB2;
        const double Q01 = -(fx*ey + ex*fy - 3.0*c*fx*fy) * invB2;
        const double Q11 = -(2.0*fy*ey + c*(1.0 - 3.0*fy*fy)) * invB2;

        const double invAB = 1.0 / (a * b);
        const double R00 = (1.0 - fx*fx - ex*ex + c*ex*fx) * invAB;
        const double R01 = (-fx*fy - ex*ey + c*ex*fy) * invAB;
        const double R10 = (-fy*fx - ey*ex + c*ey*fx) * invAB;
        const double R11 = (1.0 - fy*fy - ey*ey + c*ey*fy) * invAB;

        // K_uu, K_vv, K_uv in (u,v) coordinates.
        const double A00 = U2*gux*gux + U1*P00;
        const double A01 = U2*gux*guy + U1*P01;
        const double A10 = U2*guy*gux + U1*P01;
        const double A11 = U2*guy*guy + U1*P11;
        const double B00 = U2*gvx*gvx + U1*Q00;
        const double B01 = U2*gvx*gvy + U1*Q01;
        const double B10 = U2*gvy*gvx + U1*Q01;
        const double B11 = U2*gvy*gvy + U1*Q11;
        const double C00 = U2*gux*gvx + U1*R00;
        const double C01 = U2*gux*gvy + U1*R01;
        const double C10 = U2*guy*gvx + U1*R10;
        const double C11 = U2*guy*gvy + U1*R11;

        auto& tl = triplets_local[tid];
        const auto addEntry = [&tl, &dofMap](const int row, const int col, const double value)
        {
            const int reducedRow = dofMap[row];
            const int reducedCol = dofMap[col];
            if (reducedRow < 0 || reducedCol < 0) return;
            tl.emplace_back(reducedRow, reducedCol, value);
        };
        const auto addBlock = [&addEntry](const int rowNode, const int colNode,
                                          const double b00, const double b01,
                                          const double b10, const double b11)
        {
            addEntry(2*rowNode,     2*colNode,     b00);
            addEntry(2*rowNode,     2*colNode + 1, b01);
            addEntry(2*rowNode + 1, 2*colNode,     b10);
            addEntry(2*rowNode + 1, 2*colNode + 1, b11);
        };

        // Universal block transform of physicsDerivation.md section 5.4.
        addBlock(i, i, A00, A01, A10, A11);
        addBlock(i, j, -A00 - C00, -A01 - C01, -A10 - C10, -A11 - C11);
        addBlock(i, k, C00, C01, C10, C11);
        addBlock(j, i, -A00 - C00, -A10 - C10, -A01 - C01, -A11 - C11);
        addBlock(j, j,
                 A00 + C00 + C00 + B00,
                 A01 + C01 + C10 + B01,
                 A10 + C10 + C01 + B10,
                 A11 + C11 + C11 + B11);
        addBlock(j, k, -C00 - B00, -C01 - B10, -C10 - B01, -C11 - B11);
        addBlock(k, i, C00, C10, C01, C11);
        addBlock(k, j, -C00 - B00, -C10 - B01, -C01 - B10, -C11 - B11);
        addBlock(k, k, B00, B01, B10, B11);
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

    // Builds [K A^T; A 0] and solves it for the given K. A shifted K or the
    // identity fallback can have a different sparsity pattern, so those calls
    // force a fresh symbolic analysis.
    const auto solveAugmented = [&](const SpMat& Kmat, const bool forceAnalyze,
                                    const VecXd& violation,
                                    VecXd& dx, VecXd& dlambda) -> bool
    {
        std::vector<Eigen::Triplet<double>> augmentedTriplets;
        augmentedTriplets.reserve(Kmat.nonZeros() + 2 * A.nonZeros());

        for (int j = 0; j < Kmat.outerSize(); j++)
            for (SpMat::InnerIterator it(Kmat, j); it; ++it)
                augmentedTriplets.emplace_back(it.row(), it.col(), it.value());

        for (int j = 0; j < A.outerSize(); j++)
            for (SpMat::InnerIterator it(A, j); it; ++it)
            {
                augmentedTriplets.emplace_back(it.col(), numFreeDofs + it.row(), it.value());
                augmentedTriplets.emplace_back(numFreeDofs + it.row(), it.col(), it.value());
            }

        SpMat Kaug(total, total);
        Kaug.setFromTriplets(augmentedTriplets.begin(), augmentedTriplets.end());

        if (forceAnalyze || !patternAnalyzed)
        {
            solver.analyzePattern(Kaug);
            patternAnalyzed = true;
        }
        solver.factorize(Kaug);
        if (solver.info() != Eigen::Success) return false;

        VecXd rhs(total);
        rhs.head(numFreeDofs) = -g;
        if (numConstraints > 0) rhs.tail(numConstraints) = -violation;

        const VecXd sol = solver.solve(rhs);
        if (solver.info() != Eigen::Success) return false;

        dx = sol.head(numFreeDofs);
        dlambda = numConstraints > 0 ? sol.tail(numConstraints) : VecXd(0);
        return dx.allFinite() && dlambda.allFinite();
    };

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

        const double violationNorm = violation.norm();
        if (g.norm() < TOL && violationNorm < TOL)
        {
            std::cout << "Converged at iteration " << iter << "\n";
            break;
        }

        VecXd dx, dlambda;
        bool usedDifferentPattern = false;
        bool solved = solveAugmented(K, false, violation, dx, dlambda);
        double slope = solved ? g.dot(dx) : 0.0;

        // A moving constraint (DISPLACE_NORMAL) makes the KKT step a feasibility
        // step: the target moved, so the energy may legitimately increase. Only a
        // step with a feasible target is required to be an energy-descent step.
        const bool feasibilityStep = violationNorm > CONSTRAINT_FEASIBILITY_TOL;

        if (!solved || (!feasibilityStep && (!std::isfinite(slope) || slope >= 0.0)))
        {
            // Indefinite Hessian (or failed factorization): regularise K until the
            // step is finite (and, on a feasible target, a descent direction).
            solved = false;
            double shift = HESSIAN_SHIFT_INITIAL;
            while (shift <= HESSIAN_SHIFT_MAX)
            {
                SpMat shifted = K;
                for (int d = 0; d < numFreeDofs; d++) shifted.coeffRef(d, d) += shift;
                if (solveAugmented(shifted, true, violation, dx, dlambda))
                {
                    usedDifferentPattern = true;
                    slope = g.dot(dx);
                    if (feasibilityStep || (std::isfinite(slope) && slope < 0.0))
                    {
                        solved = true;
                        break;
                    }
                }
                shift *= HESSIAN_SHIFT_FACTOR;
            }
        }

        if (!solved)
        {
            // Last resort: projected steepest descent, still satisfying A dx = -v.
            SpMat identity(numFreeDofs, numFreeDofs);
            identity.setIdentity();
            solved = solveAugmented(identity, true, violation, dx, dlambda);
            usedDifferentPattern = true;
            slope = solved ? g.dot(dx) : 0.0;
        }

        if (usedDifferentPattern) patternAnalyzed = false;

        if (!solved || !dx.allFinite())
        {
            std::cerr << "Solid: could not compute a step at iteration " << iter
                      << "; stopping." << std::endl;
            break;
        }

        if (!feasibilityStep && !std::isfinite(slope))
        {
            std::cerr << "Solid: non-finite slope at iteration " << iter
                      << "; stopping." << std::endl;
            break;
        }
        if (!feasibilityStep
            && slope >= -1e-11 * std::max(1.0, g.norm() * dx.norm()))
        {
            // The projected gradient is zero: this is a constrained KKT point even
            // when the unprojected reduced gradient is not. Nothing left to do.
            break;
        }

        const VecXd x0 = mesh.nodePosFlat;
        VecXd delta = VecXd::Zero(mesh.nodePosFlat.size());
        for (int r = 0; r < numFreeDofs; r++)
            delta[reducedToFreeDof[r]] = dx[r];

        double alpha = 1.0;
        if (feasibilityStep)
        {
            mesh.nodePosFlat = x0 + delta;
        }
        else
        {
            // Backtracking line search on the total potential energy. Because
            // A dx = -violation, the constraint violation contracts as (1-alpha).
            const double energy0 = computeTotalPotentialEnergy();
            bool accepted = false;
            for (int ls = 0; ls < LINE_SEARCH_MAX_ITERATIONS; ls++)
            {
                mesh.nodePosFlat = x0 + alpha * delta;
                const double trialEnergy = computeTotalPotentialEnergy();
                if (std::isfinite(trialEnergy)
                    && trialEnergy <= energy0 + LINE_SEARCH_ARMIJO * alpha * slope)
                {
                    accepted = true;
                    break;
                }
                alpha *= LINE_SEARCH_REDUCTION;
                if (alpha < LINE_SEARCH_MIN_STEP) break;
            }
            if (!accepted)
            {
                mesh.nodePosFlat = x0;
                std::cerr << "Solid: line search failed at iteration " << iter
                          << " (energy " << energy0 << ", slope " << slope
                          << "); stopping." << std::endl;
                break;
            }
        }

        if (!mesh.nodePosFlat.allFinite())
        {
            mesh.nodePosFlat = x0;
            std::cerr << "Solid: non-finite step at iteration " << iter
                      << "; stopping." << std::endl;
            break;
        }

        lambda += alpha * dlambda;

        const double res_g = g.norm();
        const double res_c = violationNorm;
        const double stepNorm = alpha * dx.norm();
        if (res_g < TOL && res_c < TOL && stepNorm < STEP_TOL)
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

double Solid::computeTorsionalEnergyFromCache() const
{
    double totalEnergy = 0.0;
    for (Index t = 0; t < mesh.torsions.numTorsions; t++)
    {
        const TorsionProperties& property = mesh.torsions.torsionProperty.at(t);
        if (property.cached_a < EPSILON || property.cached_b < EPSILON) continue;
        const double deltaCosine = property.cached_cosine - property.natural_cosine;
        totalEnergy += 0.5 * property.stiffness * deltaCosine * deltaCosine;
    }
    return totalEnergy;
}

double Solid::computeTotalPotentialEnergy()
{
    mesh.cacheEdgesExtensionFromFlatVec();
    mesh.cacheTorsionStateFromFlatVec();
    return computeElasticEnergyFromExtensionCache() + computeTorsionalEnergyFromCache();
}

void Solid::assembleEnergyDerivatives()
{
    assembleGradientAndHessian();
}

void Solid::logCurrentElasticEnergy()
{
    const double energy = computeTotalPotentialEnergy();
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
