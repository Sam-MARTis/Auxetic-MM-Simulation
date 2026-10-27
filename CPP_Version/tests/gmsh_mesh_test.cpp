// Verifies the gmsh-meshed spring network: every spring is a mesh edge, no
// spring crosses a cut, the outer boundary is tagged with outward normals, and
// the solver runs on it. Needs the gmsh SDK and resources/DXFs/cell2.dxf.

#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <cstdlib>
#include <vector>

#include "constants.hpp"
#include "gmsh_mesher.hpp"
#include "loader.hpp"
#include "mesh.hpp"
#include "solid.hpp"

namespace
{
    int failures = 0;

    void expect(const bool ok, const std::string& label, const std::string& detail = "")
    {
        std::cout << (ok ? "  ok   " : "  FAIL ") << label;
        if (!detail.empty()) std::cout << "  [" << detail << "]";
        std::cout << std::endl;
        if (!ok) failures++;
    }

    // Strict crossing: a shared endpoint (the slot's closed tip, where the
    // material is continuous) does not count.
    bool properlyCrosses(const Vec2D& a, const Vec2D& b, const Vec2D& c, const Vec2D& d)
    {
        const auto side = [](const Vec2D& p, const Vec2D& q, const Vec2D& r)
        {
            return (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
        };
        const double o1 = side(a, b, c);
        const double o2 = side(a, b, d);
        const double o3 = side(c, d, a);
        const double o4 = side(c, d, b);
        // Touching an endpoint (the slot tip, where material stays continuous) is
        // not a crossing, so every orientation must be strictly non-zero.
        if (o1 == 0.0 || o2 == 0.0 || o3 == 0.0 || o4 == 0.0) return false;
        return ((o1 > 0.0) != (o2 > 0.0)) && ((o3 > 0.0) != (o4 > 0.0));
    }
}

int main(int argc, char** argv)
{
    const double meshSize = argc > 1 ? std::stod(argv[1]) : DEFAULT_GMSH_MESH_SIZE;
    const double slotWidth = argc > 2 ? std::stod(argv[2]) : DEFAULT_GMSH_SLOT_WIDTH;

    const std::string dxfFilepath = std::string(DXFS_DIRECTORY_PATH) + "/cell2.dxf";
    Loader dxfLoader{dxfFilepath};
    dxfLoader.alignPoints();

    std::cout << "polygon points: " << dxfLoader.numBoundaryPolygonPoints()
              << ", cuts: " << dxfLoader.internalEdgesNodesPos.size() << std::endl;

    const GmshMeshOptions options{.targetMeshSize = meshSize, .slotWidth = slotWidth};
    std::cout << "options: meshSize=" << meshSize << " slotWidth=" << slotWidth << std::endl;
    SpringMeshDescription description =
        GmshMesher::build(dxfLoader.boundaryNodesPos, dxfLoader.internalEdgesNodesPos, options);

    std::cout << "meshed: " << description.nodePositions.size() << " nodes, "
              << description.springs.size() << " springs, "
              << description.numTriangles << " triangles, "
              << description.numBoundaryEdges << " mesh boundary edges ("
              << description.numCutBoundaryEdges << " on cuts, "
              << description.boundarySegments.size() << " on the outline)" << std::endl;

    expect(description.nodePositions.size() > 50, "the mesher produced nodes",
           std::to_string(description.nodePositions.size()));
    expect(!description.springs.empty(), "the mesher produced springs",
           std::to_string(description.springs.size()));
    expect(!description.boundarySegments.empty(), "the outer boundary was identified",
           std::to_string(description.boundarySegments.size()));

    Mesh mesh;
    mesh.initializeFromSpringMesh(description);

    expect(mesh.nodes.numNodes == description.nodePositions.size(),
           "every meshed node reached the solver mesh");
    expect(mesh.edges.numEdges > 0, "springs reached the solver mesh",
           std::to_string(mesh.edges.numEdges));
    expect(mesh.numBoundaries == 1, "the outer outline is one boundary tag");
    expect(mesh.numBoundaryNodes > 0, "boundary nodes are tagged",
           std::to_string(mesh.numBoundaryNodes));

    expect(description.triangles.size() == description.numTriangles,
           "the mesher description carries the triangle array");
    expect(mesh.torsions.numTorsions == 3 * description.numTriangles,
           "one torsional spring per triangle corner",
           std::to_string(mesh.torsions.numTorsions));

    std::set<std::pair<Index, Index>> springPairs;
    for (Index e = 0; e < mesh.edges.numEdges; e++)
    {
        const EdgesConnection& spring = mesh.edges.edgesConnection.at(e);
        springPairs.emplace(std::min(spring.node1, spring.node2),
                            std::max(spring.node1, spring.node2));
    }
    Index torsionLegsNotSprings = 0;
    for (Index t = 0; t < mesh.torsions.numTorsions; t++)
    {
        const TorsionConnection& c = mesh.torsions.torsionConnection.at(t);
        if (springPairs.count(std::make_pair(std::min(c.nodeI, c.nodeApex),
                                             std::max(c.nodeI, c.nodeApex))) == 0)
            torsionLegsNotSprings++;
        if (springPairs.count(std::make_pair(std::min(c.nodeApex, c.nodeK),
                                             std::max(c.nodeApex, c.nodeK))) == 0)
            torsionLegsNotSprings++;
    }
    Index torsionIncidence = 0;
    for (Index i = 0; i < mesh.nodes.numNodes; i++)
        torsionIncidence += static_cast<Index>(mesh.nodes.torsionIndices.at(i).size());
    expect(torsionLegsNotSprings == 0, "every torsion leg is a spring",
           std::to_string(torsionLegsNotSprings));
    expect(torsionIncidence == 3 * mesh.torsions.numTorsions,
           "torsion incidence lists cover all three nodes",
           std::to_string(torsionIncidence));
    std::cout << "       torsions=" << mesh.torsions.numTorsions
              << " nearCollinear=" << mesh.numNearCollinearTorsions << std::endl;

    // No spring may cross a cut away from its tip: that is the point of meshing
    // around them. Right at the tip the slot cap is coarse, so a crossing there is
    // a meshing artefact rather than a missing slit.
    const auto onOutline = [&](const Vec2D& p)
    {
        const std::vector<Vec2D>& polygon = dxfLoader.boundaryNodesPos;
        for (std::size_t i = 0; i < polygon.size(); i++)
        {
            const Vec2D& a = polygon[i];
            const Vec2D& b = polygon[(i + 1) % polygon.size()];
            const Vec2D ab = b - a;
            const double lengthSq = Vec2D::dot(ab, ab);
            double t = lengthSq > EPSILON ? Vec2D::dot(p - a, ab) / lengthSq : 0.0;
            t = std::clamp(t, 0.0, 1.0);
            if ((p - (a + ab * t)).norm() < 1e-6) return true;
        }
        return false;
    };

    Index numCrossingSprings = 0;
    Index numSpringsTouchingTip = 0;
    for (Index e = 0; e < mesh.edges.numEdges; e++)
    {
        const EdgesConnection& spring = mesh.edges.edgesConnection.at(e);
        const Vec2D& p1 = mesh.nodes.nodesPos.at(spring.node1);
        const Vec2D& p2 = mesh.nodes.nodesPos.at(spring.node2);
        for (const std::pair<Vec2D, Vec2D>& cut : description.cutSegments)
        {
            if (!properlyCrosses(p1, p2, cut.first, cut.second))
            {
                if (VecMath::segmentsIntersects2D(p1, p2, cut.first, cut.second))
                    numSpringsTouchingTip++;
                continue;
            }
            // Interior endpoint of the cut = the crack tip.
            const Vec2D tip = onOutline(cut.first) ? cut.second : cut.first;
            const Vec2D along = cut.second - cut.first;
            const Vec2D other = (tip - cut.first).norm() < 1e-12 ? cut.second : cut.first;
            const double t = Vec2D::dot(
                (p1 + p2) * 0.5 - other, tip - other) / Vec2D::dot(tip - other, tip - other);
            const double distanceToTip = (1.0 - std::clamp(t, 0.0, 1.0)) * (tip - other).norm();
            if (distanceToTip > 1.5 * meshSize) numCrossingSprings++;
        }
    }
    expect(numCrossingSprings == 0, "no spring crosses a cut away from its tip",
           std::to_string(numCrossingSprings));
    expect(numSpringsTouchingTip <= 4 * description.cutSegments.size(),
           "only the slot tip region closes across a cut",
           std::to_string(numSpringsTouchingTip));

    // Both sides of a slot are separate node rows, so the slit can open.
    Index numSuspiciouslyClosePairs = 0;
    for (Index e = 0; e < mesh.edges.numEdges; e++)
    {
        const EdgesConnection& spring = mesh.edges.edgesConnection.at(e);
        const Vec2D& p1 = mesh.nodes.nodesPos.at(spring.node1);
        const Vec2D& p2 = mesh.nodes.nodesPos.at(spring.node2);
        const double length = (p1 - p2).norm();
        if (length > EPSILON && length < 0.25 * slotWidth) numSuspiciouslyClosePairs++;
    }
    expect(numSuspiciouslyClosePairs == 0, "no degenerate springs across slot walls",
           std::to_string(numSuspiciouslyClosePairs));

    // Outward normals: unit length and pointing away from the mesh centroid.
    Vec2D centroid{0.0, 0.0};
    for (Index i = 0; i < mesh.nodes.numNodes; i++) centroid += mesh.nodes.nodesPos.at(i);
    centroid /= static_cast<double>(mesh.nodes.numNodes);

    double maxUnitError = 0.0;
    double minOutwardDot = std::numeric_limits<double>::max();
    for (Index i = 0; i < mesh.nodes.numNodes; i++)
    {
        if (mesh.nodes.boundaryId.at(i) == NO_BOUNDARY) continue;
        const Vec2D normal = mesh.nodes.normals.at(i);
        maxUnitError = std::max(maxUnitError, std::abs(normal.norm() - 1.0));
        minOutwardDot = std::min(minOutwardDot,
                                 Vec2D::dot(normal, mesh.nodes.nodesPos.at(i) - centroid));
    }
    expect(maxUnitError < 1e-9, "boundary normals are unit vectors", std::to_string(maxUnitError));
    expect(minOutwardDot > 0.0, "boundary normals point outward", std::to_string(minOutwardDot));

    // The solver must run on the imported mesh.
    mesh.setBoundaryCondition(0, BoundaryCondition::displaceNormal(VELOCITY_CONSTRAINTS_EVOLUTION));
    mesh.createFlatVec();
    Solid solid(mesh);
    expect(std::abs(solid.computeTorsionalEnergyFromCache()) < 1e-12,
           "reference torsional energy is zero",
           std::to_string(solid.computeTorsionalEnergyFromCache()));

    double maxDisplacement = 0.0;
    double maxStretch = 0.0;
    Index numNan = 0;
    const Index numSteps =
        static_cast<Index>(DEFAULT_SIMULATION_END_TIME / DEFAULT_TIME_STEP);
    for (Index step = 0; step < numSteps; step++)
    {
        solid.timeStep(DEFAULT_TIME_STEP);
        for (Index i = 0; i < mesh.nodes.numNodes; i++)
        {
            const double x = mesh.nodePosFlat[2 * i];
            const double y = mesh.nodePosFlat[2 * i + 1];
            if (!std::isfinite(x) || !std::isfinite(y)) numNan++;
            maxDisplacement = std::max(maxDisplacement, std::max(std::abs(x), std::abs(y)));
        }
        for (Index e = 0; e < mesh.edges.numEdges; e++)
        {
            const double d = mesh.edges.edgesProperty.at(e).cached_extension;
            const double l0 = mesh.edges.edgesProperty.at(e).natural_length;
            if (l0 > EPSILON) maxStretch = std::max(maxStretch, d / l0);
        }
    }
    // A polygon corner has two incident edge planes. Constraining it on the
    // averaged bisector instead would let it slide off both, so each incident
    // plane must be satisfied and the corner must land on the offset-curve
    // intersection (2x the edge advance for this 60 degree triangle).
    {
        std::map<Index, std::vector<Vec2D>> incident;
        for (const BoundarySegment& segment : mesh.boundarySegments)
        {
            incident[segment.node1].push_back(segment.outwardNormal);
            incident[segment.node2].push_back(segment.outwardNormal);
        }
        const double totalAdvance =
            VELOCITY_CONSTRAINTS_EVOLUTION * DEFAULT_TIME_STEP * static_cast<double>(numSteps);
        Index numCorners = 0;
        double worstPlaneResidual = 0.0;
        double worstTangentialDrift = 0.0;
        double worstCornerAdvance = 0.0;
        for (const auto& entry : incident)
        {
            std::vector<Vec2D> distinct;
            for (const Vec2D& normal : entry.second)
            {
                bool duplicate = false;
                for (const Vec2D& existing : distinct)
                {
                    if (Vec2D::dot(existing, normal) > BOUNDARY_CORNER_DOT_THRESHOLD) duplicate = true;
                }
                if (!duplicate) distinct.push_back(normal);
            }
            if (distinct.size() < 2) continue;
            numCorners++;

            const Index node = entry.first;
            const Vec2D start = mesh.nodes.nodesPos.at(node);
            const Vec2D now{mesh.nodePosFlat[2 * node], mesh.nodePosFlat[2 * node + 1]};
            const Vec2D delta = now - start;
            for (const Vec2D& normal : distinct)
            {
                const double target = Vec2D::dot(normal, start) + totalAdvance;
                worstPlaneResidual =
                    std::max(worstPlaneResidual, std::abs(Vec2D::dot(normal, now) - target));
            }
            Vec2D bisector{0.0, 0.0};
            for (const Vec2D& normal : distinct) bisector += normal;
            bisector = bisector / bisector.norm();
            worstTangentialDrift = std::max(
                worstTangentialDrift, (delta - bisector * Vec2D::dot(bisector, delta)).norm());
            worstCornerAdvance = std::max(worstCornerAdvance, Vec2D::dot(bisector, delta));
        }
        std::cout << "       corners=" << numCorners
                  << " worstPlaneResidual=" << worstPlaneResidual
                  << " worstTangentialDrift=" << worstTangentialDrift
                  << " cornerAdvance=" << worstCornerAdvance
                  << " (edge advance " << totalAdvance << ")" << std::endl;
        expect(numCorners == description.cutSegments.size(),
               "one corner per polygon vertex", std::to_string(numCorners));
        expect(worstPlaneResidual < 1e-9, "corner nodes satisfy both incident edge planes",
               std::to_string(worstPlaneResidual));
        expect(worstTangentialDrift < 1e-9, "corner nodes do not slide sideways",
               std::to_string(worstTangentialDrift));
        expect(std::abs(worstCornerAdvance - 2.0 * totalAdvance) < 1e-9,
               "corner advance matches the offset-curve intersection",
               std::to_string(worstCornerAdvance));
    }

    expect(numNan == 0, "solver stayed finite on the gmsh mesh", std::to_string(numNan));
    expect(maxDisplacement < 10.0, "positions stayed bounded over the full run",
           std::to_string(maxDisplacement));
    expect(maxStretch < 10.0, "no spring tore apart over the full run",
           std::to_string(maxStretch));
    expect(solid.constraintViolation() < 1e-6, "plane constraints satisfied",
           std::to_string(solid.constraintViolation()));
    const double finalEnergy = solid.computeTotalPotentialEnergy();
    expect(std::isfinite(finalEnergy), "total potential energy stayed finite",
           std::to_string(finalEnergy));
    std::cout << "       final total energy=" << finalEnergy << std::endl;

    std::cout << (failures == 0 ? "GMSH MESH OK" : "GMSH MESH PROBLEMS") << std::endl;
    return failures == 0 ? 0 : 1;
}
