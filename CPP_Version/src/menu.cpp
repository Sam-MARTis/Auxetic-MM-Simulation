#include "menu.hpp"

#include "imgui-SFML.h"
#include "imgui.h"
#include "flag.hpp"


Menu::Menu(MenuData& menuDataInput): menuData(menuDataInput)
{
}

void Menu::processMenu(Renderer& renderer)
{
    StateUpdateFlags renderFlags = menuData.renderInfo.propertyUpdated;
    RendererSetupInformation& info = menuData.renderInfo.info;
    ImGui::SliderFloat("Particle Radius", &info.pRadius, 1e-3, 2e-2);
    if (ImGui::IsItemDeactivatedAfterEdit())
        renderFlags |= StateUpdateFlags::PARTICLE_RENDER_RADIUS;

    if ((renderFlags & RebuildRenderMask) != StateUpdateFlags::NONE)
    {
        renderFlags |= StateUpdateFlags::SHOULD_REBUILD_VERTICES;
        menuData.renderInfo.propertyUpdated = renderFlags;
        renderer.updateRenderInformation(menuData.renderInfo);
        menuData.renderInfo.propertyUpdated &= ~RebuildRenderMask;
    }
    ImGui::Checkbox("Render", &menuData.renderInfo.shouldRender);

    bool showNormals = renderer.showBoundaryNormalsEnabled();
    if (ImGui::Checkbox("Boundary normals", &showNormals))
    {
        renderer.setShowBoundaryNormals(showNormals);
    }
}
