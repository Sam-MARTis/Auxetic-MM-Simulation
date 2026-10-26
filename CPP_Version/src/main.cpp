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
#if defined(HAS_GMSH)
#include "gmsh_mesher.hpp"
#endif
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

    BoundaryMeshBuildInformation boundaryInfo{
        .targetDensity = TARGET_BOUNDARY_NODE_DENSITY,
        .boundary = dxfLoader.boundaryNodesPos,
    };
    boundaryInfo.boundary.printPoints();

    Mesh mesh;

    InternalEdgesMeshBuildInformation internalEdgeInfo{
        .targetDensity = TARGET_INTERNAL_EDGE_NODE_DENSITY,
        .internalEdges = dxfLoader.internalEdgesNodesPos,
    };

    const auto buildWithNodeSampling = [&]()
    {
        mesh.initializeBoundaryNodes(boundaryInfo);
        mesh.initializeInternalEdgesNodes(internalEdgeInfo);
        mesh.calculateAreaOfMeshObject(DEFAULT_NUMBER_AREA_SAMPLES);
        mesh.addInternalNodes_Uniform(DEFAULT_NUMBER_INTERNAL_NODES);
        mesh.createEdges(DEFAULT_NUM_NEIGHBOURS_FOR_EDGES, DEFAULT_DISTANCE_THRESHOLD_FOR_EDGES);
    };

    if (DEFAULT_USE_GMSH_MESHER)
    {
#if defined(HAS_GMSH)
        GmshMeshOptions gmshOptions{
            .targetMeshSize = DEFAULT_GMSH_MESH_SIZE,
            .slotWidth = DEFAULT_GMSH_SLOT_WIDTH,
            .writeMeshFile = DEFAULT_GMSH_WRITE_MESH,
            .meshFilePath = DEFAULT_GMSH_MESH_PATH,
        };
        try
        {
            mesh.initializeFromSpringMesh(
                GmshMesher::build(dxfLoader.boundaryNodesPos, dxfLoader.internalEdgesNodesPos,
                                  gmshOptions));
            std::cout << "Meshed with gmsh: " << mesh.nodes.numNodes << " nodes, "
                      << mesh.edges.numEdges << " springs" << std::endl;
        }
        catch (const std::exception& error)
        {
            std::cerr << "gmsh meshing failed (" << error.what()
                      << "); falling back to the node-sampling mesher." << std::endl;
            buildWithNodeSampling();
        }
#else
        std::cerr << "Built without gmsh; falling back to the node-sampling mesher." << std::endl;
        buildWithNodeSampling();
#endif
    }
    else
    {
        buildWithNodeSampling();
    }
    std::cout << "Edges initialized" << std::endl;

    // boundary 0 is the auto-tagged outline. Uniform expansion is the stable
    // configuration; see DeveloperReference.md section 4.5.
    mesh.setBoundaryCondition(0, BoundaryCondition::displaceNormal(VELOCITY_CONSTRAINTS_EVOLUTION));

    // CAUTION: freezing part of the outline while driving another part outward
    // collapses the springs and diverges the solve (DeveloperReference.md 4.5).
    // const Index patch = mesh.tagNodesByDomain(
    //     Cuboid(Vec2D{DEFAULT_TAG_REGION_MIN, DEFAULT_TAG_REGION_MIN},
    //            Vec2D{DEFAULT_TAG_REGION_SPLIT, DEFAULT_TAG_REGION_MAX}));
    // mesh.setBoundaryCondition(patch, BoundaryCondition::fixedNormal());

    mesh.logBoundarySummary();

    RendererSetupInformation setupInfo{
        .renderOrigin = Vec2D{DEFAULT_SCREEN_WIDTH / 2.0, DEFAULT_SCREEN_HEIGHT / 2.0},
        .nodePositions = mesh.nodes.nodesPos,
        .edgesConnections = mesh.edges.edgesConnection,
        .nodeBoundaryIds = mesh.nodes.boundaryId,
        .nodeNormals = mesh.nodes.normals
    };

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

    // The GL context must exist before the Renderer builds its sf::Texture.
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

    Renderer testRenderer{setupInfo};

    Recorder recorder;
    Menu renderMenu(menuData);

    double time = 0;
    double dt = DEFAULT_TIME_STEP;
    Solid solid(mesh);
    solid.logBoundaryConditionSummary();
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
