#include <iostream>

#include <fstream>
#include <string>
#include <cstdio>
#include <sstream>

#include<vector>
#include<tuple>

struct Vec2D
{
public:
    double x, y;
    Vec2D(double x, double y){
        this->x = x;
        this->y = y;
    }
};
    

int main()
{
    std::string command1 = "python3 readTest.py ";
    std::string dxfFilePath = "/home/odin/Academic/Research/Metamaterials/Auxetic_Mechanical_MM/Bistable/Computational_Exploration/Modules/DXF_Making/01_Trial/test_dxfs/pattern_1x1_edge-gap_0.04.dxf";
    std::string command = (command1 + dxfFilePath);
    std::cout<<command<<std::endl;
    FILE *pipe = popen(command.c_str(), "r");
    if (!pipe)
    {
        std::cerr << "Failed to start Python script";
        return 1;
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
        return 1;
    }
    unsigned int numBoundaryPoints;
    unsigned int numCutPairs;
    stream >> numBoundaryPoints;
    
    std::vector<Vec2D> boundaryPoints;
    std::vector<std::pair<Vec2D, Vec2D>> cutPairs;


    for (unsigned int i = 0; i < numBoundaryPoints; i++)
    {
        double x1, y1;
        stream >> x1 >> y1;
        boundaryPoints.emplace_back(x1, y1);
    }
    std::cout<<"Number of boundary points read: "<<numBoundaryPoints <<std::endl;
    for (unsigned  int i = 0; i < numBoundaryPoints; i++)
    {
        const Vec2D& v = boundaryPoints.at(i);
        std::cout<<v.x<<" "<<v.y<<std::endl;
    }
    
    stream>> numCutPairs;
    for(unsigned int i =0; i<numCutPairs; i++){
        double x1, y1, x2, y2;
        stream >> x1 >> y1 >> x2 >> y2;
        cutPairs.push_back(std::make_pair(Vec2D(x1, y1), Vec2D(x2, y2)));
    }

    std::cout<<"Number of cut point pairs read: "<<numCutPairs <<std::endl;
    for (unsigned  int i = 0; i < numBoundaryPoints; i++)
    {
        const std::pair<Vec2D, Vec2D>& v = cutPairs.at(i);
        std::cout<<std::get<0>(v).x<<" "<<std::get<0>(v).y<<", "<<std::get<1>(v).x<<" "<<std::get<1>(v).y<<std::endl;
    }

}