#pragma once
#include <SFML/Graphics.hpp>

#include "utils.hpp"
#include "constants.hpp"



struct RendererSetupInformation
{
    const double pRadius = DEFAULT_PARTICLE_RENDER_RADIUS;
    const sf::Color pColour = sf::Color::Blue;
    const double edgeThickness = DEFAULT_EDGE_RENDER_THICKNESS;
    const sf::Color edgeColour = sf::Color::Red;
    const double renderScale;
    const Vec2D renderOrigin;
    const std::vector<Vec2D> nodePositions;
};
enum updatedRenderProperty
{
    PARTICLE_RENDER_RADIUS,
    PARTICLE_RENDER_COLOUR,
    EDGE_RENDER_THICKNESS,
    EDGE_RENDER_COLOUR,
    RENDER_SCALE,
    RENDER_ORIGIN,
    NODE_POSITIONS,
};

struct RendererUpdateInformatiom
{
    updatedRenderProperty propertyUpdated;
    RendererSetupInformation info;
};

class Renderer
{
    double particleRadius;
    double edgeThickness;
    sf::Color particleColour = sf::Color::Blue;
    sf::Color edgeColour = sf::Color::Red;
    Vec2D renderSpacePhysicsOrigin;
    double renderScaleFromPhysics;
    std::vector<Vec2D> particlePositions;
    sf::VertexArray particleVertices;
    unsigned int numParticles;
    sf::Texture particleTexture = sf::Texture(PARTICLE_TEXTURE);
    float tw, th;
    sf::RenderStates renderState;
    public:
    Renderer(const RendererSetupInformation& setupInfo);
    void updateRenderInformation(const RendererUpdateInformatiom& updateInfo);
    Vec2D mapPhysicsToSpaceCoords(const Vec2D& physicsIn);
    Vec2D mapPhysicsToSpaceCoords(const double& physicsInX, const double& physicsInY);
    Vec2D mapPhysicsToSpaceCoords(const float& physicsInX, const float& physicsInY);

};