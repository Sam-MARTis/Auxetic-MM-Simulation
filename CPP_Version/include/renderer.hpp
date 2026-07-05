#pragma once
#include <SFML/Graphics.hpp>

#include "utils.hpp"
#include "constants.hpp"

#ifndef PARTICLE_TEXTURE_PATH
#define PARTICLE_TEXTURE_PATH "resources/textures/circle.png"
#endif

struct RendererSetupInformation
{
    double pRadius = DEFAULT_PARTICLE_RENDER_RADIUS;
    sf::Color pColour = sf::Color::Blue;
    double edgeThickness = DEFAULT_EDGE_RENDER_THICKNESS;
    sf::Color edgeColour = sf::Color::Red;
    double renderScale = DEFAULT_SCREEN_SCALE;
    Vec2D renderOrigin;
    std::vector<Vec2D>& nodePositions;
};


struct RendererUpdateData
{
    renderUpdateFlags propertyUpdated;
    RendererSetupInformation& info;
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
    sf::Texture particleTexture = sf::Texture(PARTICLE_TEXTURE_PATH);
    float tw, th;
    sf::RenderStates renderState;

    void rebuildVertexArray();
    public:
    Renderer(const RendererSetupInformation& setupInfo);
    void updateRenderInformation(const RendererUpdateData& updateInfo);
    Vec2D mapPhysicsToSpaceCoords(const Vec2D& physicsIn) const;
    Vec2D mapPhysicsToSpaceCoords(const double& physicsInX, const double& physicsInY) const;
    Vec2D mapPhysicsToSpaceCoords(const float& physicsInX, const float& physicsInY) const;
    void draw(sf::RenderWindow& window) const;
};