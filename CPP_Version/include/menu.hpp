#pragma once
#include "renderer.hpp"


struct SimulationData
{
    bool shouldTimeStep = false;
    bool shouldRender = false;
    double dt;
    double relaxationFactor;
};
struct MenuData
{
    RendererUpdateData& renderInfo;
    SimulationData& simInfo;
};



class Menu
{
private:
    MenuData& menuData;
public:
    Menu(MenuData& menuDataInput);
    void processMenu(Renderer& renderer);
};