#include "menu.hpp"

#include "imgui-SFML.h"
#include "imgui.h"


Menu::Menu(MenuData& menuDataInput): menuData(menuDataInput)
{
}

void Menu::processMenu(Renderer& renderer)
{
    RenderUpdateFlags renderFlags = RenderUpdateFlags::NONE;
    RendererSetupInformation& info = menuData.renderInfo.info;
    ImGui::SliderFloat("Particle Radius", &info.pRadius, 1e-3, 5e-1);
    if (ImGui::IsItemDeactivatedAfterEdit())
        renderFlags |= RenderUpdateFlags::PARTICLE_RENDER_RADIUS;

    if ((renderFlags & RebuildMask) != RenderUpdateFlags::NONE)
    {
        renderFlags |= RenderUpdateFlags::SHOULD_REBUILD_VERTICES;
    }
    if (renderFlags != RenderUpdateFlags::NONE)
    {
        menuData.renderInfo.propertyUpdated = renderFlags;
        renderer.updateRenderInformation(menuData.renderInfo);
    }

}
