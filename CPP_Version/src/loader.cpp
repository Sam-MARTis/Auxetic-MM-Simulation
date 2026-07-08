#include "loader.hpp"

#include <iostream>

#include "constants.hpp"
#include "utils.hpp"
#include <sstream>

Boundary::Boundary(const std::vector<Vec2D>& points)
{
    // Copy points and set the count used by isInsideBoundary
    this->boundaryPointsPos = points;
    this->numBoundaryPoints = static_cast<unsigned int>(points.size());
}

bool Boundary::isInsideBoundary(const Vec2D& testPoint) const
{
    const Vec2D rayDirection = Vec2D(1.0, 0.0);
    unsigned int numIntersections = 0;
    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        const Vec2D& boundaryPoint = boundaryPointsPos[i];
        const Vec2D segmentDirection = boundaryPointsPos[(i+1)%numBoundaryPoints] - boundaryPointsPos[i];
        // const double denom = rayDirection.x*segmentDirection.y - rayDirection.y*segmentDirection.x;
        const double denom = crossVec2D(rayDirection, segmentDirection);
        if (std::abs(denom) <= EPSILON) continue;
        const double t = crossVec2D(boundaryPoint - testPoint, segmentDirection)/denom;
        if (t < 0.0) continue;
        const double u = crossVec2D(boundaryPoint - testPoint, rayDirection)/denom;
        if (u>=0.0 && u < 1.0) numIntersections++;

    }
    if (numIntersections%2 == 0) return false;
    return true;
}

InternalEdges::InternalEdges(const std::vector<std::pair<Vec2D, Vec2D>>& points)
{
    this->internalEdgesPointsPos = points;
    this->numInternalEdges = static_cast<unsigned int>(points.size());
}


Loader::Loader(const std::string& dxfFilepath)
{
    numBoundaryPoints = 0;
    numCutPairs = 0;
    std::string command1 = "python3 ";
    std::string pyFilePath = std::string(SUBROUTINES_PATH) + "/" + "dxfParser.py ";
    // std::string dxfFilePath = "/home/odin/Academic/Research/Metamaterials/Auxetic_Mechanical_MM/Bistable/Computational_Exploration/Modules/DXF_Making/01_Trial/test_dxfs/pattern_1x1_edge-gap_0.04.dxf";
    std::string command = (command1 + pyFilePath +  dxfFilepath);
    // std::cout<<command<<std::endl;
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
    // unsigned int numBoundaryPoints;
    // unsigned int numCutPairs;
    stream >> numBoundaryPoints;
    boundaryNodesPos.reserve(numBoundaryPoints);
    // std::vector<Vec2D> boundaryPoints;
    // std::vector<std::pair<Vec2D, Vec2D>> cutPairs;


    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        double x1, y1;
        stream >> x1 >> y1;
        boundaryNodesPos.emplace_back(x1, y1);
    }


    stream>> numCutPairs;
    internalEdgesNodesPos.reserve(numBoundaryPoints);
    for(unsigned int i =0; i<numCutPairs; i++){
        double x1, y1, x2, y2;
        stream >> x1 >> y1 >> x2 >> y2;
        internalEdgesNodesPos.emplace_back(Vec2D(x1, y1), Vec2D(x2, y2));
    }
}


void Loader::alignPoints()
{
    double min_x = boundaryNodesPos.at(0).x;
    double min_y = boundaryNodesPos.at(0).y;
    double max_x = boundaryNodesPos.at(0).x;
    double max_y = boundaryNodesPos.at(0).y;
    // std::cout<<"Aligning"<<std::endl;
    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        const Vec2D& nodePos = boundaryNodesPos.at(i);
        if (nodePos.x < min_x) min_x = nodePos.x;
        if (nodePos.y < min_y) min_y = nodePos.y;
        if (nodePos.x > max_x) max_x = nodePos.x;
        if (nodePos.y > max_y) max_y = nodePos.y;
    }
    for (unsigned int i = 0; i < numCutPairs; i++)
    {
        const Vec2D& nodePos1 = std::get<0>(internalEdgesNodesPos.at(i));
        const Vec2D& nodePos2 = std::get<1>(internalEdgesNodesPos.at(i));
        if (nodePos1.x < min_x) min_x = nodePos1.x;
        if (nodePos1.y < min_y) min_y = nodePos1.y;
        if (nodePos2.x < min_x) min_x = nodePos2.x;
        if (nodePos2.y < min_y) min_y = nodePos2.y;
        if (nodePos1.x > max_x) max_x = nodePos1.x;
        if (nodePos1.y > max_y) max_y = nodePos1.y;
        if (nodePos2.x > max_x) max_x = nodePos2.x;
        if (nodePos2.y > max_y) max_y = nodePos2.y;
    }
    // std::cout<<"Align factors computed"<<std::endl;
    double x_offset = (max_x + min_x)/2.0;
    double y_offset = (max_y + min_y)/2.0;
    double x_scaling = 2.0/(max_x - min_x);
    double y_scaling = 2.0/(max_y - min_y);
    double scaling = y_scaling>x_scaling?x_scaling:y_scaling;
    const Vec2D offset = Vec2D(x_offset, y_offset);
    for (unsigned int i = 0; i < numCutPairs; i++)
    {
        Vec2D& nodePos1 = std::get<0>(internalEdgesNodesPos[i]);
        nodePos1 -= offset;
        nodePos1 *= scaling;
        Vec2D& nodePos2 = std::get<1>(internalEdgesNodesPos[i]);
        nodePos2 -= offset;
        nodePos2 *= scaling;
    }
    // std::cout<<"Internal edges aligned"<<std::endl;
    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        Vec2D& nodePos = boundaryNodesPos.at(i);
        nodePos -= offset;
        nodePos *= scaling;
    }
};
