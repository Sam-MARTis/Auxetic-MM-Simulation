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
int main()
{

    omp_set_num_threads(18);
    Eigen::setNbThreads(18);
    sf::Clock clock;
    std::cout<<"Hello World!"<<std::endl;
    // std::vector<Vec2D> boundaryNodes;
    //
    // for (double angle = 0; angle < 1.5*M_PI; angle += 0.1)
    //     boundaryNodes.push_back(Vec2D(cos(angle), sin(angle)));
    //
    // std::vector<std::pair<Vec2D, Vec2D>> internalEdges;
    // internalEdges.push_back(std::make_pair(Vec2D(0, 0), Vec2D(0.5, 0.6)));
    // internalEdges.push_back(std::make_pair(Vec2D(-0.5, 0.5), Vec2D(0.4, 0.7)));
    //
    // std::string dxfFilepath = std::string(DXFS_DIRECTORY_PATH)+ "/pattern_1x1_edge-gap_0.04.dxf";
    std::string dxfFilepath = std::string(DXFS_DIRECTORY_PATH)+ "/cell2.dxf";
    Loader dxfLoader{dxfFilepath};

    dxfLoader.alignPoints();
    Mesh mesh;

    BoundaryMeshBuildInformation boundaryInfo{
        .targetDensity = TARGET_BOUNDARY_NODE_DENSITY,
            .boundary = dxfLoader.boundaryNodesPos,
        };
    const Boundary& bound =  boundaryInfo.boundary;
    bound.printBoundaryPoints();
    // bound.testPoint({.5, .5});

    InternalEdgesMeshBuildInformation internalEdgeInfo{
    .targetDensity = TARGET_INTERNAL_EDGE_NODE_DENSITY,
    .internalEdges = dxfLoader.internalEdgesNodesPos,
    };


    mesh.initializeBoundaryNodes(boundaryInfo);
    mesh.initializeInternalEdgesNodes(internalEdgeInfo);
    double dxfArea = mesh.calculateAreaOfMeshObject(DEFAULT_NUMBER_AREA_SAMPLES);
    std::cout<<"dxfArea="<<dxfArea<<std::endl;

    // mesh.addInternalNodes_Random(DEFAULT_NUMBER_INTERNAL_NODES);

    mesh.addInternalNodes_Uniform(DEFAULT_NUMBER_INTERNAL_NODES);
    mesh.createEdges(DEFAULT_NUM_NEIGHBOURS_FOR_EDGES, DEFAULT_DISTANCE_THRESHOLD_FOR_EDGES);

    // EdgesConnection testEdge1{.node1 = 0, .node2 = 900};
    // EdgesConnection testEdge2{.node1 = 100, .node2 = 200};
    // EdgesConnection test
    // std::vector<EdgesConnection> testEdges{testEdge1, testEdge2};
    std::cout <<"Edges initialized"<<std::endl;
    RendererSetupInformation setupInfo{.nodePositions = mesh.nodes.nodesPos, .edgesConnections = mesh.edges.edgesConnection, .boundaryEdgeInfo = mesh.boundaryEdgeInfo};
    setupInfo.renderOrigin = Vec2D(DEFAULT_SCREEN_WIDTH/2.0, DEFAULT_SCREEN_HEIGHT/2.0);
    RendererUpdateData renderUpdateInfo{
        .propertyUpdated = StateUpdateFlags::NONE,
        .info = setupInfo,
        .shouldRender = true,

    };
    SimulationData simInfo{
    .dt = 0.1,
    .relaxationFactor = 1.0};

    MenuData menuData{
        .renderInfo = renderUpdateInfo,
        .simInfo = simInfo
    };

    Renderer testRenderer { setupInfo };
    sf::ContextSettings contextSettings;
    contextSettings.antiAliasingLevel = 16;
    sf::RenderWindow window = sf::RenderWindow(sf::VideoMode({static_cast<unsigned int>(round(DEFAULT_SCREEN_WIDTH)), static_cast<unsigned int>(round(DEFAULT_SCREEN_HEIGHT))}),"Testing");
    bool imguiInitSuccess = ImGui::SFML::Init(window);
    if (!imguiInitSuccess)
        std::cerr << "SFML Init failed." << std::endl;
    window.setFramerateLimit(60);

    Menu renderMenu(menuData);

    


    unsigned int iterationCount = 0;
    double time = 0;
    double dt = DEFAULT_TIME_STEP;
    Solid solid(mesh);
    
    while (window.isOpen())
    {
        PROFILE_SCOPE("Main::frame");
        window.clear(sf::Color::Black);
        while (auto event = window.pollEvent())
        {
            ImGui::SFML::ProcessEvent(window, *event);
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }
        
        ImGui::SFML::Update(window, clock.restart());
        renderMenu.processMenu(testRenderer);

        if (renderMenu.menuData.renderInfo.shouldRender)
        {
            
            testRenderer.draw(window);
        }
        for(unsigned int i=0; i<2; i++)
        {
            solid.timeStep(dt);
            time += dt;
            iterationCount++;
        }
        // std::cout <<"test particle pos "<<mesh.nodePosFlat[0] <<std::endl;
        mesh.syncVecWithFlatVecData();
        testRenderer.updateNodePositionsFromMesh(mesh);
        profiling::printLiveReport();
        ImGui::SFML::Render(window);
        window.display();
    }
    return 0;
}