#include <iostream>
#include <SFML/Graphics.hpp>

#include "imgui-SFML.h"

#include "utils.hpp"
#include "renderer.hpp"
#include "constants.hpp"
#include "solid.hpp"
#include "menu.hpp"
#include "mesh.hpp"
#include "loader.hpp"
#include <omp.h>
#include <implot.h>

int main()
{
    omp_set_num_threads(DEFAULT_NUM_THREADS);
    Eigen::setNbThreads(DEFAULT_NUM_THREADS);
    sf::Clock clock;
    std::cout << "Hello World!" << std::endl;

    std::string dxfFilepath = std::string(DXFS_DIRECTORY_PATH) + "/cell2.dxf";
    Loader dxfLoader{dxfFilepath};

    dxfLoader.alignPoints();
    Mesh mesh;

    BoundaryMeshBuildInformation boundaryInfo{
        .targetDensity = TARGET_BOUNDARY_NODE_DENSITY,
        .boundary = dxfLoader.boundaryNodesPos,
    };
    const Boundary& bound = boundaryInfo.boundary;
    bound.printBoundaryPoints();

    InternalEdgesMeshBuildInformation internalEdgeInfo{
        .targetDensity = TARGET_INTERNAL_EDGE_NODE_DENSITY,
        .internalEdges = dxfLoader.internalEdgesNodesPos,
    };

    mesh.initializeBoundaryNodes(boundaryInfo);
    mesh.initializeInternalEdgesNodes(internalEdgeInfo);
    const double dxfArea = mesh.calculateAreaOfMeshObject(DEFAULT_NUMBER_AREA_SAMPLES);
    std::cout << "dxfArea=" << dxfArea << std::endl;

    mesh.addInternalNodes_Uniform(DEFAULT_NUMBER_INTERNAL_NODES);
    mesh.createEdges(DEFAULT_NUM_NEIGHBOURS_FOR_EDGES, DEFAULT_DISTANCE_THRESHOLD_FOR_EDGES);
    std::cout << "Edges initialized" << std::endl;

    RendererSetupInformation setupInfo{
        .nodePositions = mesh.nodes.nodesPos,
        .edgesConnections = mesh.edges.edgesConnection,
        .boundaryEdgeInfo = mesh.boundaryEdgeInfo
    };
    setupInfo.renderOrigin = Vec2D{DEFAULT_SCREEN_WIDTH / 2.0, DEFAULT_SCREEN_HEIGHT / 2.0};

    RendererUpdateData renderUpdateInfo{
        .propertyUpdated = StateUpdateFlags::NONE,
        .info = setupInfo,
        .shouldRender = true,
    };
    SimulationData simInfo{
        .dt = 0.1,
        .relaxationFactor = 1.0
    };
    MenuData menuData{
        .renderInfo = renderUpdateInfo,
        .simInfo = simInfo
    };

    Renderer testRenderer{setupInfo};
    sf::ContextSettings contextSettings;
    contextSettings.antiAliasingLevel = DEFAULT_ANTIALIASING_LEVEL;
    sf::RenderWindow window = sf::RenderWindow(
        sf::VideoMode({static_cast<unsigned int>(round(DEFAULT_SCREEN_WIDTH)),
                       static_cast<unsigned int>(round(DEFAULT_SCREEN_HEIGHT))}),
        "Testing");
    if (!ImGui::SFML::Init(window))
        std::cerr << "SFML Init failed." << std::endl;
    ImPlot::CreateContext();
    window.setFramerateLimit(DEFAULT_TARGET_FRAMERATE);

    Recorder recorder;
    Menu renderMenu(menuData);

    double time = 0;
    double dt = DEFAULT_TIME_STEP;
    Solid solid(mesh);
    recorder.logEnergy(time, 0.0);

    while (window.isOpen())
    {
        if (time > DEFAULT_SIMULATION_END_TIME) break;
        PROFILE_SCOPE("Main::frame");
        window.clear(sf::Color::White);
        while (auto event = window.pollEvent())
        {
            ImGui::SFML::ProcessEvent(window, *event);
            if (event->is<sf::Event::Closed>())
                window.close();
        }

        ImGui::SFML::Update(window, clock.restart());
        ImGui::Text("Data points: %d", static_cast<int>(recorder.timeLog().size()));
        renderMenu.processMenu(testRenderer);

        if (renderMenu.menuData.renderInfo.shouldRender)
            testRenderer.draw(window);

        recorder.recordFrame(window);

        for (unsigned int i = 0; i < DEFAULT_SIM_STEPS_PER_FRAME; i++)
        {
            solid.timeStep(dt);
            time += dt;
            recorder.logEnergy(time, solid.computeElasticEnergyFromExtensionCache());
        }

        ImGui::Begin("Energy vs Time");
        ImPlot::SetNextAxesToFit();
        if (ImPlot::BeginPlot("Energy over Time"))
        {
            ImPlot::SetupAxes("time", "energy");
            ImPlot::PlotLine("Energy",
                             recorder.timeLog().data(),
                             recorder.energyLog().data(),
                             static_cast<int>(recorder.timeLog().size()));
            ImPlot::EndPlot();
        }
        ImGui::End();

        mesh.syncVecWithFlatVecData();
        testRenderer.updateNodePositionsFromMesh(mesh);
        ImGui::SFML::Render(window);
        window.display();
    }

    recorder.finalize();
    ImPlot::DestroyContext();

    return 0;
}
