#pragma once
#include "renderer.hpp"


struct SimulationData
{
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
public:
    MenuData& menuData;
    Menu(MenuData& menuDataInput);
    void processMenu(Renderer& renderer);
};