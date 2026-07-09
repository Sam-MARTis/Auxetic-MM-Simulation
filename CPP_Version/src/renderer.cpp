
#include "renderer.hpp"

#include <iostream>

#include "utils.hpp"

void Renderer::rebuildParticleVertexArray()
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

    const std::vector<sf::Color> colours = {
        sf::Color::Green,
        sf::Color::Red,
        sf::Color::White,
        sf::Color::Yellow,
        sf::Color::Magenta,
        sf::Color::Cyan,
        // sf::Color::Black
    };
    const unsigned int numBoundaries = boundaryEdgeInfos.size();
    for (unsigned int i=0; i<numBoundaries; i++)
    {
        const BoundaryEdgeInformation& boundaryInfo = boundaryEdgeInfos.at(i);
        const unsigned int startIndex = std::get<0>(boundaryInfo.boundsIndices);
        // std::cout<<startIndex<<std::endl;
        const unsigned int endIndex = std::get<1>(boundaryInfo.boundsIndices);
        const sf::Color thisEdgeColour = colours.at(i);
        for (unsigned int j=startIndex; j<endIndex; j++)
        {
            const Vec2D screenPos = mapPhysicsToSpaceCoords(particlePositions[j]);
            const float x = static_cast<float>(screenPos.x);
            const float y = static_cast<float>(screenPos.y);
            // Triangle 1
            particleVertices[6*j + 0] = {{x - r, y - r}, thisEdgeColour, {0.f, 0.f}};
            particleVertices[6*j + 1] = {{x - r, y + r}, thisEdgeColour, {0.f, th}};
            particleVertices[6*j + 2] = {{x + r, y - r}, thisEdgeColour, {tw, 0.f}};

            // Triangle 2
            particleVertices[6*j + 3] = {{x - r, y + r}, thisEdgeColour, {0.f, th}};
            particleVertices[6*j + 4] = {{x + r, y + r}, thisEdgeColour, {tw, th}};
            particleVertices[6*j + 5] = {{x + r, y - r}, thisEdgeColour, {tw, 0.f}};
        }
    }

}


void Renderer::rebuildEdgeVertexArray()
{
    assert(numEdges == edgeConnections.size());
    for (unsigned int i=0; i<numEdges; i++)
    {
        const EdgesConnection& edge = edgeConnections.at(i);
        const Vec2D& node1Pos = particlePositions.at(edge.node1);
        const Vec2D& node2Pos = particlePositions.at(edge.node2);
        const Vec2D screenPos1 = mapPhysicsToSpaceCoords(node1Pos);
        const Vec2D screenPos2 = mapPhysicsToSpaceCoords(node2Pos);
        const Vec2D edgeVector = screenPos2 - screenPos1;
        const Vec2D Normal = {edgeVector.y, -edgeVector.x}.normalized();
        const float halfThickness = 0.5f*edgeThickness;
        const Vec2D offset = Normal*halfThickness;
        const Vec2D v1 = screenPos1 + offset;
        const Vec2D v2 = screenPos1 - offset;
        const Vec2D v3 = screenPos2 + offset;
        const Vec2D v4 = screenPos2 - offset;
        // Triangle 1
        edgeVertices[6*i + 0] = {{static_cast<float>(v1.x), static_cast<float>(v1.y)}, edgeColour};
        edgeVertices[6*i + 1] = {{static_cast<float>(v2.x), static_cast<float>(v2.y)}, edgeColour};
        edgeVertices[6*i + 2] = {{static_cast<float>(v3.x), static_cast<float>(v3.y)}, edgeColour};
        // Triangle 2
        edgeVertices[6*i + 3] = {{static_cast<float>(v2.x), static_cast<float>(v2.y)}, edgeColour};
        edgeVertices[6*i + 4] = {{static_cast<float>(v4.x), static_cast<float>(v4.y)}, edgeColour};
        edgeVertices[6*i + 5] = {{static_cast<float>(v3.x), static_cast<float>(v3.y)}, edgeColour};
    }
}

Renderer::Renderer(const RendererSetupInformation& setupInfo)
: boundaryEdgeInfos(setupInfo.boundaryEdgeInfo)
{
    particleRadius = setupInfo.pRadius;
    edgeThickness = setupInfo.edgeThickness;
    particleColour = setupInfo.pColour;
    edgeThickness = setupInfo.edgeThickness;
    edgeColour = setupInfo.edgeColour;
    renderScaleFromPhysics = setupInfo.renderScale;
    renderSpacePhysicsOrigin = setupInfo.renderOrigin;

    numParticles = setupInfo.nodePositions.size();
    numEdges = setupInfo.edgesConnections.size();
    particleVertices = sf::VertexArray(sf::PrimitiveType::Triangles, 6*numParticles);
    edgeVertices = sf::VertexArray(sf::PrimitiveType::Triangles, 6*numEdges);
    const sf::Vector2u texSize = particleTexture.getSize();
    tw = static_cast<float>(texSize.x);
    th = static_cast<float>(texSize.y);

    renderState.texture = &particleTexture;

    const std::vector<Vec2D>& pos = setupInfo.nodePositions;
    const std::vector<EdgesConnection>& edgeCon = setupInfo.edgesConnections;
    particlePositions.reserve(numParticles);
    particlePositions.insert(particlePositions.end(), pos.begin(), pos.end());

    edgeConnections.reserve(numEdges);
    edgeConnections.insert(edgeConnections.end(), edgeCon.begin(), edgeCon.end());
    rebuildParticleVertexArray();
    rebuildEdgeVertexArray();

}

void Renderer::updateRenderInformation(const RendererUpdateData& updateInfo)
{

    const RendererSetupInformation& info = updateInfo.info;
    StateUpdateFlags flags = updateInfo.propertyUpdated;
    if ((flags & StateUpdateFlags::NODE_POSITIONS) != StateUpdateFlags::NONE)
    {
        particlePositions.clear();
        particlePositions.insert(
            particlePositions.end(),
            info.nodePositions.begin(),
            info.nodePositions.end()
        );
    }

    if ((flags & StateUpdateFlags::RENDER_SCALE) != StateUpdateFlags::NONE)
    {
        renderScaleFromPhysics = info.renderScale;
    }

    if ((flags & StateUpdateFlags::RENDER_ORIGIN) != StateUpdateFlags::NONE)
    {
        renderSpacePhysicsOrigin = info.renderOrigin;
    }

    if ((flags & StateUpdateFlags::PARTICLE_RENDER_RADIUS) != StateUpdateFlags::NONE)
    {
        particleRadius = info.pRadius;
    }

    if ((flags & StateUpdateFlags::PARTICLE_RENDER_COLOUR) != StateUpdateFlags::NONE)
    {
        particleColour = info.pColour;
    }

    if ((flags & StateUpdateFlags::EDGE_RENDER_THICKNESS) != StateUpdateFlags::NONE)
    {
        edgeThickness = info.edgeThickness;
    }

    if ((flags & StateUpdateFlags::EDGE_RENDER_COLOUR) != StateUpdateFlags::NONE)
    {
        edgeColour = info.edgeColour;
    }

    if ((flags & StateUpdateFlags::SHOULD_REBUILD_VERTICES) != StateUpdateFlags::NONE)
    {
        rebuildParticleVertexArray();

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
    window.draw(edgeVertices, renderState);
}
