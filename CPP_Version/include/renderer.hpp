#pragma once
#include <SFML/Graphics.hpp>

#include <filesystem>
#include <string>
#include <vector>

#include "utils.hpp"
#include "constants.hpp"
#include "flag.hpp"
#include "mesh.hpp"

#ifndef PARTICLE_TEXTURE_PATH
#define PARTICLE_TEXTURE_PATH "resources/textures/circle.png"
#endif

struct RendererSetupInformation
{
    float pRadius = DEFAULT_PARTICLE_RENDER_RADIUS;
    sf::Color pColour = sf::Color::Blue;
    float edgeThickness = DEFAULT_EDGE_RENDER_THICKNESS;
    sf::Color edgeColour = sf::Color::Red;
    float renderScale = DEFAULT_SCREEN_SCALE;
    Vec2D renderOrigin;
    std::vector<Vec2D>& nodePositions;
    std::vector<EdgesConnection>& edgesConnections;

    // Per-node boundary tag (Mesh::nodes::boundaryId); NO_BOUNDARY means interior.
    std::vector<Index>& nodeBoundaryIds;
    // Per-node outward unit normal, drawn as a short segment from tagged nodes.
    std::vector<Vec2D>& nodeNormals;
    // Off by default: the per-node ticks are useful for checking that "outward"
    // really is outward, but they clutter the mesh view.
    bool showBoundaryNormals = false;
    float normalLength = DEFAULT_BOUNDARY_NORMAL_RENDER_LENGTH;
    sf::Color normalColour = sf::Color::Black;
};


struct RendererUpdateData
{
    StateUpdateFlags propertyUpdated;
    RendererSetupInformation& info;
    bool shouldRender;
};

class Renderer
{
    float particleRadius;
    float edgeThickness;
    sf::Color particleColour = sf::Color::Blue;
    sf::Color edgeColour = sf::Color::Red;
    Vec2D renderSpacePhysicsOrigin;
    float renderScaleFromPhysics;
    std::vector<Vec2D> particlePositions;
    std::vector<EdgesConnection> edgeConnections;
    std::vector<Index> nodeBoundaryIds;
    std::vector<Vec2D> nodeNormals;
    sf::VertexArray particleVertices;
    sf::VertexArray edgeVertices;
    sf::VertexArray normalVertices;
    float normalLength;
    sf::Color normalColour = sf::Color::Black;
    bool showBoundaryNormals = false;
    unsigned int numParticles;
    unsigned int numEdges;
    sf::Texture particleTexture = sf::Texture(PARTICLE_TEXTURE_PATH);
    // A 1x1 opaque white texture used for geometry that must not sample the
    // particle sprite. It exists so that NO draw call ever passes a null
    // texture: see the note on Renderer::draw.
    sf::Texture opaqueTexture;
    float tw, th;
    sf::RenderStates renderState;
    sf::RenderStates opaqueRenderState;

    void rebuildParticleVertexArray();
    void rebuildEdgeVertexArray();
    void rebuildNormalVertexArray();
    public:
    Renderer(const RendererSetupInformation& setupInfo);
    void updateRenderInformation(const RendererUpdateData& updateInfo);
    Vec2D mapPhysicsToSpaceCoords(const Vec2D& physicsIn) const;
    Vec2D mapPhysicsToSpaceCoords(const double& physicsInX, const double& physicsInY) const;
    Vec2D mapPhysicsToSpaceCoords(const float& physicsInX, const float& physicsInY) const;
    void draw(sf::RenderWindow& window) const;
    void updateNodePositionsFromMesh(Mesh& mesh);
    void setShowBoundaryNormals(bool show);
    [[nodiscard]] bool showBoundaryNormalsEnabled() const;
};

struct RecorderSetupInformation
{
    bool recordStructureEvolution = true;
    bool snapshotEnergyLog = true;
    unsigned int videoFrameRate = DEFAULT_RECORDING_FRAME_RATE;
    std::filesystem::path framesDirectory = DEFAULT_RECORDING_FRAMES_DIR;
    std::filesystem::path videosDirectory = OUTPUT_VIDEOS_DIR;
    std::filesystem::path imagesDirectory = OUTPUT_IMAGES_DIR;
    std::filesystem::path subroutinesDirectory = SUBROUTINES_PATH;
};

// Captures simulation frames, energy history, the energy log snapshot and the
// assembled structure-evolution video. Everything it writes lives under the
// directories configured in RecorderSetupInformation.
class Recorder
{
    RecorderSetupInformation setup;
    std::vector<double> times;
    std::vector<double> energies;
    unsigned int recordedFrameCount = 0;

    void writeEnergyLog() const;
    void createEnergySnapshot() const;
    void encodeVideo() const;
public:
    explicit Recorder(const RecorderSetupInformation& setupInfo = {});
    void recordFrame(sf::RenderWindow& window);
    void logEnergy(double time, double energy);
    void finalize();
    [[nodiscard]] const std::vector<double>& timeLog() const;
    [[nodiscard]] const std::vector<double>& energyLog() const;
};
