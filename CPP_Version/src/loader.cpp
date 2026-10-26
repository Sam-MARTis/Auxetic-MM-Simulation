#include "loader.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>

#include "constants.hpp"
#include "utils.hpp"

namespace
{
    // Chains unordered boundary segments into one closed polygon, or returns an
    // empty vector when they do not form a single closed loop.
    std::vector<Vec2D> assemblePolygon(const std::vector<std::pair<Vec2D, Vec2D>>& lines)
    {
        if (lines.empty()) return {};

        double minX = lines[0].first[0], maxX = minX;
        double minY = lines[0].first[1], maxY = minY;
        for (const std::pair<Vec2D, Vec2D>& line : lines)
        {
            for (const Vec2D& p : {line.first, line.second})
            {
                minX = std::min(minX, p[0]);
                maxX = std::max(maxX, p[0]);
                minY = std::min(minY, p[1]);
                maxY = std::max(maxY, p[1]);
            }
        }
        const double scale = std::max(maxX - minX, maxY - minY);
        const double tolerance = 1e-9 * std::max(scale, 1.0);

        std::vector<bool> used(lines.size(), false);
        std::vector<Vec2D> polygon;
        polygon.reserve(lines.size());
        polygon.push_back(lines[0].first);
        Vec2D current = lines[0].second;
        used[0] = true;

        for (std::size_t step = 1; step < lines.size(); step++)
        {
            polygon.push_back(current);
            bool advanced = false;
            for (std::size_t i = 0; i < lines.size() && !advanced; i++)
            {
                if (used[i]) continue;
                if ((lines[i].first - current).norm() <= tolerance)
                {
                    current = lines[i].second;
                    used[i] = true;
                    advanced = true;
                }
                else if ((lines[i].second - current).norm() <= tolerance)
                {
                    current = lines[i].first;
                    used[i] = true;
                    advanced = true;
                }
            }
            if (!advanced) return {};
        }
        if ((current - polygon.front()).norm() > tolerance) return {};
        return polygon;
    }
}

InternalEdgesDescription::InternalEdgesDescription(const std::vector<std::pair<Vec2D, Vec2D>>& points)
{
    this->internalEdgesPointsPos = points;
    this->numInternalEdges = static_cast<unsigned int>(points.size());
}

Loader::Loader(const std::string& dxfFilepath)
{
    numBoundaryLines = 0;
    numCutPairs = 0;
    std::string command1 = "python3 ";
    std::string pyFilePath = std::string(SUBROUTINES_PATH) + "/" + "dxfParser.py ";
    std::string command = (command1 + pyFilePath + dxfFilepath);
    FILE *pipe = popen(command.c_str(), "r");
    if (!pipe)
    {
        std::cerr << "Failed to start Python script";
        exit(1);
    }
    std::stringstream stream;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe))
    {
        stream << buffer;
    }

    int status = pclose(pipe);
    if (status != 0)
    {
        std::cerr << "Python script failed.\n";
    }

    stream >> numBoundaryLines;
    std::vector<std::pair<Vec2D, Vec2D>> boundaryLines;
    boundaryLines.reserve(numBoundaryLines);
    for (unsigned int i = 0; i < numBoundaryLines; i++)
    {
        double x1, y1, x2, y2;
        stream >> x1 >> y1 >> x2 >> y2;
        boundaryLines.emplace_back(Vec2D{x1, y1}, Vec2D{x2, y2});
    }

    stream >> numCutPairs;
    internalEdgesNodesPos.reserve(numCutPairs);
    for (unsigned int i = 0; i < numCutPairs; i++)
    {
        double x1, y1, x2, y2;
        stream >> x1 >> y1 >> x2 >> y2;
        internalEdgesNodesPos.emplace_back(Vec2D{x1, y1}, Vec2D{x2, y2});
    }

    boundaryNodesPos = assemblePolygon(boundaryLines);
    if (boundaryNodesPos.empty())
    {
        // Fall back to the segment start points, which is what the DXF produces
        // when the lines are already ordered end to end.
        std::cerr << "Loader: boundary lines do not form a closed loop, "
                     "falling back to their start points." << std::endl;
        boundaryNodesPos.reserve(boundaryLines.size());
        for (const std::pair<Vec2D, Vec2D>& line : boundaryLines) boundaryNodesPos.push_back(line.first);
    }
}

void Loader::alignPoints()
{
    const std::size_t numBoundaryPoints = boundaryNodesPos.size();
    double min_x = boundaryNodesPos.at(0)[0];
    double min_y = boundaryNodesPos.at(0)[1];
    double max_x = boundaryNodesPos.at(0)[0];
    double max_y = boundaryNodesPos.at(0)[1];
    for (std::size_t i = 0; i < numBoundaryPoints; i++)
    {
        const Vec2D& nodePos = boundaryNodesPos.at(i);
        if (nodePos[0] < min_x) min_x = nodePos[0];
        if (nodePos[1] < min_y) min_y = nodePos[1];
        if (nodePos[0] > max_x) max_x = nodePos[0];
        if (nodePos[1] > max_y) max_y = nodePos[1];
    }
    for (unsigned int i = 0; i < numCutPairs; i++)
    {
        const Vec2D& nodePos1 = std::get<0>(internalEdgesNodesPos.at(i));
        const Vec2D& nodePos2 = std::get<1>(internalEdgesNodesPos.at(i));
        if (nodePos1[0] < min_x) min_x = nodePos1[0];
        if (nodePos1[1] < min_y) min_y = nodePos1[1];
        if (nodePos2[0] < min_x) min_x = nodePos2[0];
        if (nodePos2[1] < min_y) min_y = nodePos2[1];
        if (nodePos1[0] > max_x) max_x = nodePos1[0];
        if (nodePos1[1] > max_y) max_y = nodePos1[1];
        if (nodePos2[0] > max_x) max_x = nodePos2[0];
        if (nodePos2[1] > max_y) max_y = nodePos2[1];
    }
    double x_offset = (max_x + min_x)/2.0;
    double y_offset = (max_y + min_y)/2.0;
    double x_scaling = 2.0/(max_x - min_x);
    double y_scaling = 2.0/(max_y - min_y);
    double scaling = y_scaling>x_scaling?x_scaling:y_scaling;
    const Vec2D offset = Vec2D{x_offset, y_offset};
    for (std::size_t i = 0; i < numBoundaryPoints; i++)
    {
        Vec2D& nodePos = boundaryNodesPos.at(i);
        nodePos -= offset;
        nodePos *= scaling;
    }
    for (unsigned int i = 0; i < numCutPairs; i++)
    {
        Vec2D& nodePos1 = std::get<0>(internalEdgesNodesPos[i]);
        nodePos1 -= offset;
        nodePos1 *= scaling;
        Vec2D& nodePos2 = std::get<1>(internalEdgesNodesPos[i]);
        nodePos2 -= offset;
        nodePos2 *= scaling;
    }
}
