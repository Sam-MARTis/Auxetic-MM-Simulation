
#include "renderer.hpp"
#include "utils.hpp"

void Renderer::rebuildVertexArray()
{
    const float r = renderScaleFromPhysics*static_cast<float>(particleRadius);
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
    const sf::Vector2u texSize = particleTexture.getSize();
    tw = static_cast<float>(texSize.x);
    th = static_cast<float>(texSize.y);

    renderState.texture = &particleTexture;

    const std::vector<Vec2D>& pos = setupInfo.nodePositions;
    particlePositions.reserve(numParticles);
    particlePositions.insert(particlePositions.end(), pos.begin(), pos.end());
    rebuildVertexArray();
}

void Renderer::updateRenderInformation(const RendererUpdateData& updateInfo)
{
    const RendererSetupInformation& info = updateInfo.info;
    RenderUpdateFlags flags = updateInfo.propertyUpdated;
    if ((flags & RenderUpdateFlags::NODE_POSITIONS) != RenderUpdateFlags::NONE)
    {
        particlePositions.clear();
        particlePositions.insert(
            particlePositions.end(),
            info.nodePositions.begin(),
            info.nodePositions.end()
        );
    }

    if ((flags & RenderUpdateFlags::RENDER_SCALE) != RenderUpdateFlags::NONE)
    {
        renderScaleFromPhysics = info.renderScale;
    }

    if ((flags & RenderUpdateFlags::RENDER_ORIGIN) != RenderUpdateFlags::NONE)
    {
        renderSpacePhysicsOrigin = info.renderOrigin;
    }

    if ((flags & RenderUpdateFlags::PARTICLE_RENDER_RADIUS) != RenderUpdateFlags::NONE)
    {
        particleRadius = info.pRadius;
    }

    if ((flags & RenderUpdateFlags::PARTICLE_RENDER_COLOUR) != RenderUpdateFlags::NONE)
    {
        particleColour = info.pColour;
    }

    if ((flags & RenderUpdateFlags::EDGE_RENDER_THICKNESS) != RenderUpdateFlags::NONE)
    {
        edgeThickness = info.edgeThickness;
    }

    if ((flags & RenderUpdateFlags::EDGE_RENDER_COLOUR) != RenderUpdateFlags::NONE)
    {
        edgeColour = info.edgeColour;
    }

    if ((flags & RenderUpdateFlags::SHOULD_REBUILD_VERTICES) != RenderUpdateFlags::NONE)
    {
        rebuildVertexArray();
    }
}

Vec2D Renderer::mapPhysicsToSpaceCoords(const Vec2D& physicsIn) const
{
    return {renderSpacePhysicsOrigin.x + renderScaleFromPhysics*physicsIn.x, renderSpacePhysicsOrigin.y - renderScaleFromPhysics* physicsIn.y};
}

Vec2D Renderer::mapPhysicsToSpaceCoords(const double& physicsInX, const double& physicsInY) const
{
    return  {renderSpacePhysicsOrigin.x + renderScaleFromPhysics*physicsInX, renderSpacePhysicsOrigin.y - renderScaleFromPhysics* physicsInY};
}
Vec2D Renderer::mapPhysicsToSpaceCoords(const float& physicsInX, const float& physicsInY) const
{
    return  {renderSpacePhysicsOrigin.x + renderScaleFromPhysics*physicsInX, renderSpacePhysicsOrigin.y - renderScaleFromPhysics* physicsInY};
}

void Renderer::draw(sf::RenderWindow& window) const
{
    window.draw(particleVertices, renderState);
}
