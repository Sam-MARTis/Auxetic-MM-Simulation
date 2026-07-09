
#pragma once

#include "mesh.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>


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
    private:
    std::vector<BoundaryConstraint> constraints;

    // Remember to call flat vec initialization of mesh when solid made
    void intializeNormalDirichletBCs();
    public:
    Solid(Mesh& solidMesh);
    Mesh& mesh;
    void assembleConstraintMatrix();

    void buildConstraintRHS();
    void assembleGradientAndHessian();
    void solve();
    void timeStep(const double dt);

};
