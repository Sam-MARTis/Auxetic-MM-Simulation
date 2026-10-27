#include "gmsh_mesher.hpp"

#include <gmsh.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace
{
    struct EdgeRecord
    {
        int triangleCount = 0;
        std::size_t firstTriangle = 0;
    };

    std::pair<Index, Index> orderedPair(const Index a, const Index b)
    {
        return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
    }

    double distanceToSegment(const Vec2D& p, const Vec2D& a, const Vec2D& b)
    {
        const Vec2D ab = b - a;
        const double lengthSq = Vec2D::dot(ab, ab);
        if (lengthSq <= EPSILON) return (p - a).norm();
        const double t = std::clamp(Vec2D::dot(p - a, ab) / lengthSq, 0.0, 1.0);
        return (p - (a + ab * t)).norm();
    }

    bool onPolygon(const Vec2D& p, const std::vector<Vec2D>& polygon, const double tolerance)
    {
        const std::size_t n = polygon.size();
        for (std::size_t i = 0; i < n; i++)
        {
            if (distanceToSegment(p, polygon[i], polygon[(i + 1) % n]) <= tolerance) return true;
        }
        return false;
    }

    // Adds a closed polygon of straight edges and returns its plane surface tag.
    int addPolygonSurface(const std::vector<Vec2D>& polygon, const double meshSize)
    {
        const std::size_t numPoints = polygon.size();
        std::vector<int> pointTags(numPoints);
        std::vector<int> lineTags;
        lineTags.reserve(numPoints);
        for (std::size_t i = 0; i < numPoints; i++)
        {
            pointTags[i] = gmsh::model::occ::addPoint(polygon[i][0], polygon[i][1], 0.0, meshSize);
        }
        for (std::size_t i = 0; i < numPoints; i++)
        {
            lineTags.push_back(
                gmsh::model::occ::addLine(pointTags[i], pointTags[(i + 1) % numPoints]));
        }
        const int loopTag = gmsh::model::occ::addCurveLoop(lineTags);
        return gmsh::model::occ::addPlaneSurface({loopTag});
    }

    // The four corners of the slot covering [from, to], widened by the slot width
    // and extended past `from` so the boolean cleanly opens the outline.
    std::vector<Vec2D> slotPolygon(const Vec2D& from, const Vec2D& to,
                                   const double slotWidth, const double extension)
    {
        Vec2D direction = to - from;
        const double length = direction.norm();
        if (length <= EPSILON) return {};
        direction /= length;
        const Vec2D normal{-direction[1], direction[0]};
        const double halfWidth = 0.5 * std::max(slotWidth, 1e-4);
        const Vec2D start = from - direction * extension;
        return {
            start + normal * halfWidth,
            to + normal * halfWidth,
            to - normal * halfWidth,
            start - normal * halfWidth,
        };
    }
}

SpringMeshDescription GmshMesher::build(const std::vector<Vec2D>& boundaryPolygon,
                                        const std::vector<std::pair<Vec2D, Vec2D>>& cuts,
                                        const GmshMeshOptions& options)
{
    if (boundaryPolygon.size() < 3)
    {
        throw std::runtime_error("GmshMesher: the boundary polygon needs at least 3 points");
    }

    // An exact slit would model the cut better, but gmsh cannot make the 2D mesh
    // conform to a curve ending on the outline (embed hangs, fragment does not
    // conform), so cuts are removed as slots.
    double slotWidth = options.slotWidth;
    if (!(slotWidth > 0.0))
    {
        std::cerr << "GmshMesher: slotWidth must be positive; using the mesh size." << std::endl;
        slotWidth = options.targetMeshSize;
    }

    SpringMeshDescription description;
    description.cutSegments = cuts;

    const bool ownsGmsh = gmsh::isInitialized() == 0;
    if (ownsGmsh)
    {
        gmsh::initialize(0, nullptr, false);
    }
    gmsh::option::setNumber("General.Terminal", options.verbosity);

    try
    {
        gmsh::model::add("spring_domain");

        const int surfaceTag = addPolygonSurface(boundaryPolygon, options.targetMeshSize);

        if (!cuts.empty())
        {
            gmsh::model::occ::synchronize();
            gmsh::vectorpair tools;
            for (const std::pair<Vec2D, Vec2D>& cut : cuts)
            {
                const std::vector<Vec2D> slot = slotPolygon(cut.first, cut.second, slotWidth,
                                                            options.slotExtension);
                if (slot.empty()) continue;
                tools.emplace_back(2, addPolygonSurface(slot, options.targetMeshSize));
            }

            gmsh::vectorpair result;
            std::vector<gmsh::vectorpair> resultMap;
            gmsh::model::occ::cut({{2, surfaceTag}}, tools, result, resultMap, -1, true, true);
            if (result.empty())
            {
                throw std::runtime_error("GmshMesher: the boolean cut produced no surface");
            }
        }

        gmsh::option::setNumber("Mesh.MeshSizeMax", options.targetMeshSize);
        gmsh::option::setNumber("Mesh.MeshSizeMin", 0.5 * options.targetMeshSize);
        gmsh::model::occ::synchronize();
        gmsh::model::mesh::generate(2);

        if (options.writeMeshFile)
        {
            gmsh::write(options.meshFilePath);
        }

        // ---- read the mesh back ------------------------------------------------
        std::vector<std::size_t> nodeTags;
        std::vector<double> coordinates;
        std::vector<double> parametricCoordinates;
        gmsh::model::mesh::getNodes(nodeTags, coordinates, parametricCoordinates);

        std::map<std::size_t, Index> localOf;
        description.nodePositions.reserve(nodeTags.size());
        for (std::size_t i = 0; i < nodeTags.size(); i++)
        {
            localOf[nodeTags[i]] = static_cast<Index>(i);
            description.nodePositions.emplace_back(Vec2D{coordinates[3 * i], coordinates[3 * i + 1]});
        }

        std::vector<std::size_t> elementTags;
        std::vector<std::size_t> elementNodes;
        gmsh::model::mesh::getElementsByType(2, elementTags, elementNodes);
        const std::size_t numTriangles = elementNodes.size() / 3;
        std::vector<std::array<Index, 3>> triangles(numTriangles);

        // Every triangle edge is a spring; an edge used by a single triangle is a
        // boundary of the 2D mesh.
        std::map<std::pair<Index, Index>, EdgeRecord> edges;
        for (std::size_t t = 0; t < numTriangles; t++)
        {
            for (int c = 0; c < 3; c++) triangles[t][c] = localOf.at(elementNodes[3 * t + c]);
            for (int c = 0; c < 3; c++)
            {
                EdgeRecord& record = edges[orderedPair(triangles[t][c], triangles[t][(c + 1) % 3])];
                record.triangleCount++;
                if (record.triangleCount == 1) record.firstTriangle = t;
            }
        }

        const auto oppositeCorner = [&](const std::pair<Index, Index>& edge, const std::size_t triangle)
        {
            for (int c = 0; c < 3; c++)
            {
                const Index corner = triangles[triangle][c];
                if (corner != edge.first && corner != edge.second) return corner;
            }
            return edge.first;
        };

        description.numTriangles = static_cast<Index>(numTriangles);
        description.springs.reserve(edges.size());
        for (const auto& entry : edges)
        {
            const std::pair<Index, Index>& edge = entry.first;
            const EdgeRecord& record = entry.second;
            description.springs.push_back(EdgesConnection{edge.first, edge.second});

            if (record.triangleCount != 1) continue;
            description.numBoundaryEdges++;

            const Vec2D& p1 = description.nodePositions.at(edge.first);
            const Vec2D& p2 = description.nodePositions.at(edge.second);
            const Vec2D midpoint = (p1 + p2) * 0.5;

            // An edge lying on the original outline is the outer boundary;
            // everything else (slot walls and caps) is a cut.
            if (!onPolygon(midpoint, boundaryPolygon, 1e-7))
            {
                description.numCutBoundaryEdges++;
                continue;
            }

            Vec2D direction = p2 - p1;
            Vec2D normal{-direction[1], direction[0]};
            const double magnitude = normal.norm();
            if (magnitude <= EPSILON) continue;
            normal /= magnitude;
            const Vec2D& inside = description.nodePositions.at(
                oppositeCorner(edge, record.firstTriangle));
            if (Vec2D::dot(normal, midpoint - inside) < 0.0) normal *= -1.0;

            description.boundarySegments.push_back(
                BoundarySegment{edge.first, edge.second, NO_BOUNDARY, normal});
        }

        // Torsional springs are built from the triangle corners downstream.
        description.triangles = triangles;
    }
    catch (...)
    {
        if (ownsGmsh) gmsh::finalize();
        throw;
    }

    if (ownsGmsh) gmsh::finalize();
    return description;
}
