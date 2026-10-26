// Headless diagnostics for the real DXF pipeline: runs the full simulation under
// several boundary conditions and reports constraint violation, |gradient|,
// energy, stretch and position bounds per step. Build with
// -DBUILD_DXF_DIAGNOSTICS=ON; see DeveloperReference.md section 6.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "loader.hpp"
#include "mesh.hpp"
#include "solid.hpp"
#include "constants.hpp"

namespace
{
    int failures = 0;

    void report(const bool ok, const std::string& label, const std::string& detail = "")
    {
        std::cout << (ok ? "  ok   " : "  FAIL ") << label;
        if (!detail.empty()) std::cout << "  [" << detail << "]";
        std::cout << std::endl;
        if (!ok) failures++;
    }

    void warn(const std::string& label, const std::string& detail)
    {
        std::cout << "  WARN " << label << "  [" << detail << "]" << std::endl;
    }

    Vec2D flatPosition(const Mesh& mesh, const Index nodeIdx)
    {
        return Vec2D{mesh.nodePosFlat[2 * nodeIdx], mesh.nodePosFlat[2 * nodeIdx + 1]};
    }

    bool isFinite(const double value)
    {
        return std::isfinite(value);
    }

    // A cut DXF may legitimately yield several components.
    std::vector<Index> connectedComponentSizes(const Mesh& mesh)
    {
        std::vector<int> component(mesh.nodes.numNodes, -1);
        std::vector<Index> sizes;
        std::vector<Index> stack;
        for (Index seed = 0; seed < mesh.nodes.numNodes; seed++)
        {
            if (component[seed] != -1) continue;
            const int id = static_cast<int>(sizes.size());
            Index count = 0;
            component[seed] = id;
            stack.push_back(seed);
            while (!stack.empty())
            {
                const Index current = stack.back();
                stack.pop_back();
                count++;
                for (const Index e : mesh.nodes.edgeIndices.at(current))
                {
                    const EdgesConnection& c = mesh.edges.edgesConnection.at(e);
                    const Index other = (c.node1 == current) ? c.node2 : c.node1;
                    if (component[other] != -1) continue;
                    component[other] = id;
                    stack.push_back(other);
                }
            }
            sizes.push_back(count);
        }
        return sizes;
    }

    struct ScenarioResult
    {
        double worstViolation = 0.0;
        double worstAbsPos = 0.0;
        double worstStretch = 0.0;
        double worstShrink = 0.0;
        double peakEnergy = 0.0;
        bool anyNaN = false;
        Index divergenceStep = NO_BOUNDARY;   // first step whose energy jumps by >10x
    };

    Mesh buildMesh()
    {
        // addInternalNodes_Uniform uses rand(); seed so scenarios share a mesh.
        std::srand(41);
        const std::string dxfFilepath = std::string(DXFS_DIRECTORY_PATH) + "/cell2.dxf";
        Loader dxfLoader{dxfFilepath};
        dxfLoader.alignPoints();

        BoundaryMeshBuildInformation boundaryInfo{
            .targetDensity = TARGET_BOUNDARY_NODE_DENSITY,
            .boundary = dxfLoader.boundaryNodesPos,
        };
        InternalEdgesMeshBuildInformation internalEdgeInfo{
            .targetDensity = TARGET_INTERNAL_EDGE_NODE_DENSITY,
            .internalEdges = dxfLoader.internalEdgesNodesPos,
        };

        Mesh mesh;
        mesh.initializeBoundaryNodes(boundaryInfo);
        mesh.initializeInternalEdgesNodes(internalEdgeInfo);
        mesh.calculateAreaOfMeshObject(DEFAULT_NUMBER_AREA_SAMPLES);
        mesh.addInternalNodes_Uniform(DEFAULT_NUMBER_INTERNAL_NODES);
        mesh.createEdges(DEFAULT_NUM_NEIGHBOURS_FOR_EDGES, DEFAULT_DISTANCE_THRESHOLD_FOR_EDGES);
        return mesh;
    }

    ScenarioResult runScenario(const std::string& label, const bool usePatch,
                              const BoundaryCondition& patchCondition, const bool verbose)
    {
        std::cout << "=== scenario: " << label << " ===" << std::endl;
        Mesh mesh = buildMesh();

        mesh.setBoundaryCondition(0, BoundaryCondition::displaceNormal(VELOCITY_CONSTRAINTS_EVOLUTION));
        if (usePatch)
        {
            const Index patchBoundary = mesh.tagNodesByDomain(
                Cuboid(Vec2D{DEFAULT_TAG_REGION_MIN, DEFAULT_TAG_REGION_MIN},
                       Vec2D{DEFAULT_TAG_REGION_SPLIT, DEFAULT_TAG_REGION_MAX}));
            if (!mesh.boundaryNodes(patchBoundary).empty())
                mesh.setBoundaryCondition(patchBoundary, patchCondition);
        }

        const std::vector<Index> componentSizes = connectedComponentSizes(mesh);
        std::cout << "  nodes=" << mesh.nodes.numNodes << " edges=" << mesh.edges.numEdges
                  << " taggedNodes=" << mesh.numBoundaryNodes
                  << " boundaries=" << mesh.numBoundaries
                  << " graphComponents=" << componentSizes.size() << " (sizes:";
        for (const Index size : componentSizes) std::cout << " " << size;
        std::cout << ")" << std::endl;

        Vec2D centroid{0.0, 0.0};
        for (Index i = 0; i < mesh.nodes.numNodes; i++) centroid += mesh.nodes.nodesPos.at(i);
        centroid /= static_cast<double>(mesh.nodes.numNodes);

        // Adjacent normals must not oppose, or the chord inference is bogus.
        {
            double worstAdjacentDot = 1.0;
            Index worstA = NO_BOUNDARY, worstB = NO_BOUNDARY;
            for (const BoundarySegment& segment : mesh.boundarySegments)
            {
                const Vec2D na = mesh.nodes.normals.at(segment.node1);
                const Vec2D nb = mesh.nodes.normals.at(segment.node2);
                const double dot = Vec2D::dot(na, nb);
                if (dot < worstAdjacentDot)
                {
                    worstAdjacentDot = dot;
                    worstA = segment.node1;
                    worstB = segment.node2;
                }
            }
            std::cout << "  segments=" << mesh.boundarySegments.size()
                      << " worstAdjacentNormalDot=" << worstAdjacentDot
                      << " (nodes " << worstA << "," << worstB << ")" << std::endl;
            if (worstAdjacentDot < 0.0)
                warn("adjacent boundary normals oppose", "inferred chords may be wrong");
            // Per-boundary normal spread and plane offsets.
            for (Index b = 0; b < mesh.numBoundaries; b++)
            {
                if (mesh.boundaryNodes(b).empty()) continue;
                double minDot = 1.0, maxDot = -1.0;
                for (const Index i : mesh.boundaryNodes(b))
                {
                    const Vec2D n = mesh.nodes.normals.at(i);
                    const double dotCentroid =
                        Vec2D::dot(n, mesh.nodes.nodesPos.at(i) - centroid);
                    minDot = std::min(minDot, dotCentroid);
                    maxDot = std::max(maxDot, dotCentroid);
                }
                std::cout << "  boundary " << b << ": nodes=" << mesh.boundaryNodes(b).size()
                          << " outwardDotRange=[" << minDot << "," << maxDot
                          << "] cond=" << boundaryConditionTypeName(mesh.boundaryCondition(b).type)
                          << " v=" << mesh.boundaryCondition(b).velocity << std::endl;
            }
        }

        mesh.createFlatVec();
        Solid solid(mesh);

        const double dt = DEFAULT_TIME_STEP;
        const Index numSteps =
            (static_cast<Index>(DEFAULT_SIMULATION_END_TIME / (dt * DEFAULT_SIM_STEPS_PER_FRAME)) + 2)
            * DEFAULT_SIM_STEPS_PER_FRAME;

        if (verbose)
        {
            std::cout << std::setw(6) << "step" << std::setw(13) << "energy" << std::setw(13) << "maxDispl"
                      << std::setw(13) << "maxAbsPos" << std::setw(13) << "maxStretch"
                      << std::setw(13) << "minStretch" << std::setw(13) << "planeViol"
                      << std::setw(13) << "|grad|" << std::setw(6) << "nan" << std::endl;
        }

        ScenarioResult result;
        std::vector<Vec2D> previous;
        previous.reserve(mesh.nodes.numNodes);
        for (Index i = 0; i < mesh.nodes.numNodes; i++) previous.push_back(flatPosition(mesh, i));

        double previousEnergy = 0.0;
        for (Index step = 0; step < numSteps; step++)
        {
            solid.timeStep(dt);
            mesh.cacheEdgesExtensionFromFlatVec();

            Index numNan = 0;
            double maxDispl = 0.0, maxAbsPos = 0.0, maxStretch = 0.0, minStretch = 1e300;
            for (Index i = 0; i < mesh.nodes.numNodes; i++)
            {
                const Vec2D now = flatPosition(mesh, i);
                if (!isFinite(now[0]) || !isFinite(now[1])) numNan++;
                maxAbsPos = std::max(maxAbsPos, std::max(std::abs(now[0]), std::abs(now[1])));
                const Vec2D delta = now - previous[i];
                if (isFinite(delta[0]) && isFinite(delta[1])) maxDispl = std::max(maxDispl, delta.norm());
                previous[i] = now;
            }
            for (Index e = 0; e < mesh.edges.numEdges; e++)
            {
                const double d = mesh.edges.edgesProperty.at(e).cached_extension;
                const double l0 = mesh.edges.edgesProperty.at(e).natural_length;
                if (l0 <= EPSILON) continue;
                maxStretch = std::max(maxStretch, d / l0);
                minStretch = std::min(minStretch, d / l0);
            }

            const double energy = solid.computeElasticEnergyFromExtensionCache();
            const double violation = solid.constraintViolation();

            result.worstViolation = std::max(result.worstViolation, violation);
            result.worstAbsPos = std::max(result.worstAbsPos, maxAbsPos);
            result.worstStretch = std::max(result.worstStretch, maxStretch);
            result.worstShrink = std::max(result.worstShrink, 1.0 / std::max(minStretch, 1e-12));
            result.peakEnergy = std::max(result.peakEnergy, energy);
            result.anyNaN = result.anyNaN || numNan > 0;
            if (step > 0 && energy > 10.0 * std::max(previousEnergy, 1e-9) &&
                result.divergenceStep == NO_BOUNDARY)
            {
                result.divergenceStep = step;
            }
            previousEnergy = energy;

            if (verbose && (step < 12 || step % 5 == 0 || numNan > 0 || step == result.divergenceStep))
            {
                std::cout << std::setw(6) << step << std::setw(13) << energy
                          << std::setw(13) << maxDispl << std::setw(13) << maxAbsPos
                          << std::setw(13) << maxStretch << std::setw(13) << minStretch
                          << std::setw(13) << violation << std::setw(13) << solid.gradientNorm()
                          << std::setw(6) << numNan << std::endl;
            }
        }

        std::cout << "  peakEnergy=" << result.peakEnergy
                  << " worstViolation=" << result.worstViolation
                  << " worstStretch=" << result.worstStretch
                  << " worstShrink=" << result.worstShrink
                  << " worstAbsPos=" << result.worstAbsPos << std::endl;
        // Correctness, as opposed to taste about the resulting shape.
        report(result.worstViolation < 1e-6,
               "plane constraints satisfied every step", std::to_string(result.worstViolation));
        report(!result.anyNaN, "no NaN positions");
        // Stability is data; correctness is what must hold.
        if (!(result.worstAbsPos < 10.0))
            warn("positions left the screen", std::to_string(result.worstAbsPos));
        if (!(result.worstStretch < 10.0))
            warn("lattice grossly distorted", std::to_string(result.worstStretch) + "x stretch");
        if (result.divergenceStep != NO_BOUNDARY)
            warn("energy diverged",
                 "first >10x jump at step " + std::to_string(result.divergenceStep));
        return result;
    }
}

int main()
{
    std::cout << "=== DXF pipeline diagnostics ===" << std::endl;
    const double v = VELOCITY_CONSTRAINTS_EVOLUTION;

    struct NamedScenario
    {
        std::string label;
        bool usePatch;
        BoundaryCondition patch;
    };
    const std::vector<NamedScenario> scenarios = {
        {"A: uniform DISPLACE_NORMAL (pre-refactor behaviour)", false, BoundaryCondition::free()},
        {"B: patch FIXED_NORMAL (the previous demo)", true, BoundaryCondition::fixedNormal()},
        {"C: patch FIXED (fully pinned)", true, BoundaryCondition::fixed()},
        {"D: patch DISPLACE_NORMAL inward (-v)", true, BoundaryCondition::displaceNormal(-v)},
        {"E: patch DISPLACE_NORMAL outward (v/2)", true, BoundaryCondition::displaceNormal(0.5 * v)},
    };

    std::vector<ScenarioResult> results;
    for (const NamedScenario& scenario : scenarios)
    {
        results.push_back(runScenario(scenario.label, scenario.usePatch, scenario.patch, false));
    }

    std::cout << "=== comparison ===" << std::endl;
    std::cout << std::left << std::setw(52) << "scenario" << std::setw(14) << "peakEnergy"
              << std::setw(14) << "worstStretch" << std::setw(14) << "worstAbsPos"
              << std::setw(12) << "divergedAt" << std::endl;
    for (std::size_t i = 0; i < scenarios.size(); i++)
    {
        std::cout << std::left << std::setw(52) << scenarios[i].label
                  << std::setw(14) << results[i].peakEnergy
                  << std::setw(14) << results[i].worstStretch
                  << std::setw(14) << results[i].worstAbsPos
                  << std::setw(12) << (results[i].divergenceStep == NO_BOUNDARY
                                       ? std::string("no") : std::to_string(results[i].divergenceStep))
                  << std::endl;
    }

    Index numStable = 0;
    for (const ScenarioResult& result : results)
        if (result.divergenceStep == NO_BOUNDARY && result.worstAbsPos < 10.0 &&
            result.worstStretch < 10.0)
            numStable++;

    std::cout << (failures == 0 ? "DIAGNOSTICS OK" : "DIAGNOSTICS FOUND PROBLEMS")
              << " -- " << numStable << "/" << results.size()
              << " boundary-condition scenarios stayed stable" << std::endl;
    return failures == 0 ? 0 : 1;
}
