
#pragma once

#include "mesh.hpp"
#include <Eigen/Dense>
#include "profiling.hpp"
#include <Eigen/Sparse>

#include <string>
#include <unordered_map>
#include <vector>
#include <Eigen/UmfPackSupport>


struct BoundaryConstraint
{
    Vec2D nVec;
    double c;
    unsigned int boundaryId;
};

class Solid
{
    SpMat A;
    Vec c;
    SpMat K;
    Vec g;
    SpMat K_aug;
    std::vector<Eigen::Triplet<double>> triplets;
private:
    std::vector<BoundaryConstraint> constraints;

    // Remember to call flat vec initialization of mesh when solid made
    // void intializeNormalDirichletBCs();
    void assembleConstraintMatrix();
    void buildConstraintRHS();
    void assembleGradientAndHessian();
    void solve();
    bool patternAnalyzed = false;
public:
    // Eigen::SparseLU<SpMat> solver;
    Eigen::UmfPackLU<SpMat> solver;
    Solid(Mesh& solidMesh);
    Mesh& mesh;

    void timeStep(const double dt);

};
