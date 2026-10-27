#pragma once

#include "vec.hpp"
#include "mesh.hpp"
#include <Eigen/Dense>
#include "profiling.hpp"
#include <Eigen/Sparse>

#include <string>
#include <unordered_map>
#include <vector>
#include <Eigen/UmfPackSupport>

// Static energy minimiser. Conditions reduce to either DOF elimination or a
// Lagrange plane constraint n . x = c(t); see DeveloperReference.md section 3.
class Solid
{
    SpMat A;              // M x DOF_free : normal-plane constraint rows
    VecXd c;              // M           : current plane offsets
    SpMat K;
    VecXd g;
    std::vector<Eigen::Triplet<double>> triplets;

    // DOF bookkeeping. Sizes are the full 2*numNodes / reduced DOF count.
    std::vector<int> freeDofToReduced;    // global DOF -> reduced index, -1 when fixed
    std::vector<int> reducedToFreeDof;    // reduced index -> global DOF

    // One plane per node, because a spring boundary can curve.
    struct NormalConstraint
    {
        Index boundaryId = NO_BOUNDARY;
        Index nodeIndex = NO_BOUNDARY;
        Vec2D normal{0.0, 0.0};
        double velocity = 0.0;   // signed rate of c along +n (DISPLACE_NORMAL only)
    };
    std::vector<NormalConstraint> normalConstraints;
    Index numFixedNodes = 0;

    void classifyDegreesOfFreedom();
    void buildNormalConstraints();
    void assembleConstraintMatrix();
    void buildConstraintRHS();
    void assembleGradientAndHessian();
    void applyBoundaryVelocity(double dt);
    void solve();
    [[nodiscard]] VecXd gatherReducedPositions() const;

    bool patternAnalyzed = false;
    bool nearCollinearWarned = false;
public:
    Eigen::UmfPackLU<SpMat> solver;
    Solid(Mesh& solidMesh);
    std::vector<double> EnergyHistory;
    Mesh& mesh;
    double computeElasticEnergyFromExtensionCache();
    // Reads the torsion cache; call after Mesh::cacheTorsionStateFromFlatVec().
    [[nodiscard]] double computeTorsionalEnergyFromCache() const;
    // Refreshes both derived caches and returns extension + torsional energy.
    // This is the line-search merit function.
    double computeTotalPotentialEnergy();
    // Refreshes g and K at the current nodePosFlat (diagnostics and tests).
    void assembleEnergyDerivatives();
    void timeStep(const double dt);
    void logCurrentElasticEnergy();
    // Required after changing Mesh::boundaryConditions post-construction.
    void refreshBoundaryConditions();
    [[nodiscard]] Index numConstrainedNodes() const { return static_cast<Index>(normalConstraints.size()); }
    [[nodiscard]] Index numFixedNodesCount() const { return numFixedNodes; }
    // Reduced gradient / Hessian of the total potential energy.
    [[nodiscard]] const VecXd& gradient() const { return g; }
    [[nodiscard]] const SpMat& hessian() const { return K; }
    // max |n . x - c| over plane rows (0 if none); ~1e-16 after a good solve.
    [[nodiscard]] double constraintViolation() const;
    [[nodiscard]] double gradientNorm() const;
    void logBoundaryConditionSummary() const;
};
