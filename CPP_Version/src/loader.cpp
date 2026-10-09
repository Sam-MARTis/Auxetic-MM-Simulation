#include "loader.hpp"

#include <iostream>

#include "constants.hpp"
#include "utils.hpp"
#include <sstream>
#include <cassert>
InternalEdgesDescription::InternalEdgesDescription(const std::vector<std::pair<Vec2D, Vec2D>>& points)
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
        boundaryNodesPos.emplace_back(Vec2D{x1, y1});
    }


    stream>> numCutPairs;
    internalEdgesNodesPos.reserve(numBoundaryPoints);
    for(unsigned int i =0; i<numCutPairs; i++){
        double x1, y1, x2, y2;
        stream >> x1 >> y1 >> x2 >> y2;
        internalEdgesNodesPos.emplace_back(Vec2D{x1, y1}, Vec2D{x2, y2});
    }
}


void Loader::alignPoints()
{
    double min_x = boundaryNodesPos.at(0)[0];
    double min_y = boundaryNodesPos.at(0)[1];
    double max_x = boundaryNodesPos.at(0)[0];
    double max_y = boundaryNodesPos.at(0)[1];
    // std::cout<<"Aligning"<<std::endl;
    for (unsigned int i = 0; i < numBoundaryPoints; i++)
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
    // std::cout<<"Align factors computed"<<std::endl;
    double x_offset = (max_x + min_x)/2.0;
    double y_offset = (max_y + min_y)/2.0;
    double x_scaling = 2.0/(max_x - min_x);
    double y_scaling = 2.0/(max_y - min_y);
    double scaling = y_scaling>x_scaling?x_scaling:y_scaling;
    const Vec2D offset = Vec2D{x_offset, y_offset};
    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        Vec2D& nodePos = boundaryNodesPos.at(i);
        nodePos -= offset;
        nodePos *= scaling;
    }
    scaling *= DEFAULT_INTERNAL_EDGE_EXPANSIVE_FACTOR;
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
};
