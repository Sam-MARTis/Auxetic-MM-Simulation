
#include "renderer.hpp"
#include "solid.hpp"

#include <SFML/OpenGL.hpp>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "utils.hpp"

void Renderer::rebuildParticleVertexArray()
{
    PROFILE_SCOPE("Renderer::rebuildParticleVertexArray");

    const float r = renderScaleFromPhysics*static_cast<float>(particleRadius);
    for (unsigned int i = 0; i < numParticles; i++)
    {
        const Vec2D screenPos = mapPhysicsToSpaceCoords(particlePositions[i]);
        const float x = static_cast<float>(screenPos[0]);
        const float y = static_cast<float>(screenPos[1]);
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
            const float x = static_cast<float>(screenPos[0]);
            const float y = static_cast<float>(screenPos[1]);
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
    PROFILE_SCOPE("Renderer::rebuildEdgeVertexArray");
    assert(numEdges == edgeConnections.size());
    for (unsigned int i=0; i<numEdges; i++)
    {
        const EdgesConnection& edge = edgeConnections.at(i);
        const Vec2D& node1Pos = particlePositions.at(edge.node1);
        const Vec2D& node2Pos = particlePositions.at(edge.node2);
        const Vec2D edgeVec = node2Pos - node1Pos;
        const Vec2D normal = Vec2D{edgeVec[1], -edgeVec[0]}.normalized();
        const float halfThickness = 0.5f*edgeThickness;
        const Vec2D offset = normal*halfThickness;
        const Vec2D r1 = node1Pos + offset;
        const Vec2D r2 = node1Pos - offset;
        const Vec2D r3 = node2Pos + offset;
        const Vec2D r4 = node2Pos - offset;
        const Vec2D v1 = mapPhysicsToSpaceCoords(r1);
        const Vec2D v2 = mapPhysicsToSpaceCoords(r2);
        const Vec2D v3 = mapPhysicsToSpaceCoords(r3);
        const Vec2D v4 = mapPhysicsToSpaceCoords(r4);
        // const Vec2D screenPos1 = mapPhysicsToSpaceCoords(node1Pos);
        // const Vec2D screenPos2 = mapPhysicsToSpaceCoords(node2Pos);
        // const Vec2D edgeVector = screenPos2 - screenPos1;
        // const Vec2D Normal = Vec2D{edgeVector[1], -edgeVector[0]}.normalized();
        // const float halfThickness = 0.5f*edgeThickness;
        // const Vec2D offset = Normal*halfThickness;
        // const Vec2D v1 = screenPos1 + offset;
        // const Vec2D v2 = screenPos1 - offset;
        // const Vec2D v3 = screenPos2 + offset;
        // const Vec2D v4 = screenPos2 - offset;
        // Triangle 1
        edgeVertices[6*i + 0] = {{static_cast<float>(v1[0]), static_cast<float>(v1[1])}, edgeColour, {0.4f*tw, 0.4f*th}};
        edgeVertices[6*i + 1] = {{static_cast<float>(v2[0]), static_cast<float>(v2[1])}, edgeColour, {0.4f*tw, 0.6f*th}};
        edgeVertices[6*i + 2] = {{static_cast<float>(v3[0]), static_cast<float>(v3[1])}, edgeColour, {0.6f*tw, 0.4f*th}};
        // Triangle 2
        edgeVertices[6*i + 3] = {{static_cast<float>(v2[0]), static_cast<float>(v2[1])}, edgeColour, {0.4f*tw, 0.6f*th}};
        edgeVertices[6*i + 4] = {{static_cast<float>(v4[0]), static_cast<float>(v4[1])}, edgeColour, {0.6f*tw, 0.6f*th}};
        edgeVertices[6*i + 5] = {{static_cast<float>(v3[0]), static_cast<float>(v3[1])}, edgeColour, {0.6f*tw, 0.4f*th}};
    }
}

Renderer::Renderer(const RendererSetupInformation& setupInfo)
: boundaryEdgeInfos(setupInfo.boundaryEdgeInfo)
{
    PROFILE_SCOPE("Renderer::Renderer");
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
    PROFILE_SCOPE("Renderer::updateRenderInformation");

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
    return Vec2D{renderSpacePhysicsOrigin[0] + renderScaleFromPhysics*physicsIn[0], renderSpacePhysicsOrigin[1] - renderScaleFromPhysics* physicsIn[1]};
}

Vec2D Renderer::mapPhysicsToSpaceCoords(const double& physicsInX, const double& physicsInY) const
{
    return  {renderSpacePhysicsOrigin[0] + renderScaleFromPhysics*physicsInX, renderSpacePhysicsOrigin[1] - renderScaleFromPhysics* physicsInY};
}
Vec2D Renderer::mapPhysicsToSpaceCoords(const float& physicsInX, const float& physicsInY) const
{
    return  {renderSpacePhysicsOrigin[0] + renderScaleFromPhysics*physicsInX, renderSpacePhysicsOrigin[1] - renderScaleFromPhysics* physicsInY};
}

void Renderer::draw(sf::RenderWindow& window) const
{
    PROFILE_SCOPE("Renderer::draw");
    window.draw(edgeVertices, renderState);
    window.draw(particleVertices, renderState);
    // std::cout<<numEdges<<std::endl;
}

void Renderer::updateNodePositionsFromMesh(Mesh& mesh)
{
    PROFILE_SCOPE("Renderer::updateNodePositionsFromMesh");
    particlePositions = mesh.nodes.nodesPos;
    rebuildParticleVertexArray();
    rebuildEdgeVertexArray();
}

Recorder::Recorder(const RecorderSetupInformation& setupInfo)
: setup(setupInfo)
{
    PROFILE_SCOPE("Recorder::Recorder");
    if (setup.recordStructureEvolution)
    {
        std::filesystem::remove_all(setup.framesDirectory);
        std::filesystem::create_directories(setup.framesDirectory);
        std::filesystem::create_directories(setup.videosDirectory);
    }
    if (setup.recordStructureEvolution || setup.snapshotEnergyLog)
    {
        std::filesystem::create_directories(setup.imagesDirectory);
    }
}

void Recorder::recordFrame(sf::RenderWindow& window)
{
    PROFILE_SCOPE("Recorder::recordFrame");
    if (!setup.recordStructureEvolution) return;

    const sf::Vector2u windowSize = window.getSize();
    std::vector<std::uint8_t> pixels(windowSize.x * windowSize.y * 4);
    glReadPixels(0, 0,
                 static_cast<GLsizei>(windowSize.x),
                 static_cast<GLsizei>(windowSize.y),
                 GL_RGBA,
                 GL_UNSIGNED_BYTE,
                 pixels.data());

    sf::Image frameImage({windowSize.x, windowSize.y}, pixels.data());
    frameImage.flipVertically();

    const auto framePath = setup.framesDirectory /
        (std::string(DEFAULT_RECORDING_FRAME_PREFIX) + std::to_string(recordedFrameCount) + ".png");
    if (!frameImage.saveToFile(framePath))
        std::cerr << "Failed to save recording frame: " << framePath << std::endl;
    else
        recordedFrameCount++;
}

void Recorder::logEnergy(const double time, const double energy)
{
    times.push_back(time);
    energies.push_back(energy);
}

void Recorder::writeEnergyLog() const
{
    const auto energyLogCsvPath = setup.imagesDirectory / DEFAULT_ENERGY_LOG_FILENAME;
    std::ofstream energyLogCsvFile(energyLogCsvPath);
    if (!energyLogCsvFile)
    {
        std::cerr << "Failed to open energy log CSV for writing: " << energyLogCsvPath << std::endl;
        return;
    }

    energyLogCsvFile << "time,energy\n";
    for (std::size_t i = 0; i < times.size() && i < energies.size(); ++i)
        energyLogCsvFile << times[i] << ',' << energies[i] << '\n';
}

void Recorder::createEnergySnapshot() const
{
    const auto snapshotPath = setup.imagesDirectory / DEFAULT_ENERGY_SNAPSHOT_FILENAME;
    const auto energyLogCsvPath = setup.imagesDirectory / DEFAULT_ENERGY_LOG_FILENAME;
    const std::string pythonCommand =
        std::string(DEFAULT_PYTHON_INTERPRETER) + " \"" +
        (setup.subroutinesDirectory / ENERGY_PLOT_SCRIPT_FILENAME).string() +
        "\" \"" + energyLogCsvPath.string() + "\" \"" + snapshotPath.string() + "\"";
    std::system(pythonCommand.c_str());
}

void Recorder::encodeVideo() const
{
    const auto outputVideoPath = setup.videosDirectory / DEFAULT_STRUCTURE_VIDEO_FILENAME;
    const std::string framePattern =
        (setup.framesDirectory / (std::string(DEFAULT_RECORDING_FRAME_PREFIX) + "%d.png")).string();
    const std::string ffmpegCommand =
        "ffmpeg -y -framerate " + std::to_string(setup.videoFrameRate) + " -i \"" + framePattern +
        "\" -c:v libx264 -pix_fmt yuv420p \"" + outputVideoPath.string() + "\"";
    std::system(ffmpegCommand.c_str());
}

void Recorder::finalize()
{
    PROFILE_SCOPE("Recorder::finalize");
    if (setup.snapshotEnergyLog)
    {
        writeEnergyLog();
        createEnergySnapshot();
    }
    if (setup.recordStructureEvolution)
    {
        encodeVideo();
    }
}

const std::vector<double>& Recorder::timeLog() const
{
    return times;
}

const std::vector<double>& Recorder::energyLog() const
{
    return energies;
}
