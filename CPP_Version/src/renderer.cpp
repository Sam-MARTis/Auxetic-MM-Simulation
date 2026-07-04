
#include "renderer.hpp"


void Renderer::rebuildVertexArray()
{
    const float r = static_cast<float>(particleRadius);
    for (unsigned int i = 0; i < numParticles; i++)
    {
        const Vec2D screenPos = mapPhysicsToSpaceCoords(particlePositions[i]);
        const float x = static_cast<float>(screenPos.x);
        const float y = static_cast<float>(screenPos.y);
        // Triangle 1
        particleVertices[6*i + 0] = {{x - r, y - r}, particleColour, {0.f, 0.f}};
        particleVertices[6*i + 1] = {{x - r, y + r}, particleColour, {0.f, th}};
        particleVertices[6*i + 2] = {{x + r, y - r}, particleColour, {tw, 0.f}};

        // Triangle 2
        particleVertices[6*i + 3] = {{x - r, y + r}, particleColour, {0.f, th}};
        particleVertices[6*i + 4] = {{x + r, y + r}, particleColour, {tw, th}};
        particleVertices[6*i + 5] = {{x + r, y - r}, particleColour, {tw, 0.f}};
    }
}

Renderer::Renderer(const RendererSetupInformation& setupInfo)
{
    particleRadius = setupInfo.pRadius;
    edgeThickness = setupInfo.edgeThickness;
    particleColour = setupInfo.pColour;
    edgeThickness = setupInfo.edgeThickness;
    edgeColour = setupInfo.edgeColour;
    renderScaleFromPhysics = setupInfo.renderScale;
    renderSpacePhysicsOrigin = setupInfo.renderOrigin;
    numParticles = setupInfo.nodePositions.size();
    particleVertices = sf::VertexArray(sf::PrimitiveType::Triangles, 6*numParticles);
    sf::Vector2u texSize = particleTexture.getSize();
    tw = static_cast<float>(texSize.x);
    th = static_cast<float>(texSize.y);

    renderState.texture = &particleTexture;

    const std::vector<Vec2D>& pos = setupInfo.nodePositions;
    particlePositions.reserve(numParticles);
    particlePositions.insert(particlePositions.end(), pos.begin(), pos.end());
    rebuildVertexArray();
}

void Renderer::updateRenderInformation(const RendererUpdateInformatiom& updateInfo, bool shouldRebuildVertices)
{
    const RendererSetupInformation& info = updateInfo.info;
    switch (updateInfo.propertyUpdated)
    {
        case NODE_POSITIONS:
        particlePositions.clear();
        particlePositions.insert(particlePositions.end(), info.nodePositions.begin(), info.nodePositions.end());
        break;
        case RENDER_SCALE:
        renderScaleFromPhysics = info.renderScale;
        break;
        case RENDER_ORIGIN:
        renderSpacePhysicsOrigin = info.renderOrigin;
        break;
        case PARTICLE_RENDER_RADIUS:
        particleRadius = info.pRadius;
        break;
        case PARTICLE_RENDER_COLOUR:
        particleColour = info.pColour;
        break;
        case  EDGE_RENDER_THICKNESS:
        edgeThickness = info.edgeThickness;
        break;
        case EDGE_RENDER_COLOUR:
        edgeColour = info.edgeColour;
        break;

    }
    if (shouldRebuildVertices) rebuildVertexArray();
}