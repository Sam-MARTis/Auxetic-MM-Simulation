#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>

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
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

int main()
{

    omp_set_num_threads(1);
    Eigen::setNbThreads(1);
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
    ImPlot::CreateContext();
    window.setFramerateLimit(60);





    bool should_record_structure_evolution = true;
    bool shoould_snapshot_energy_log = true;
    
    
    
    const std::filesystem::path recordingFramesDir{"recorded_structure_frames"};
    const std::filesystem::path outputVideosDir{OUTPUT_VIDEOS_DIR};
    if (should_record_structure_evolution)
    {
        std::filesystem::remove_all(recordingFramesDir);
        std::filesystem::create_directories(recordingFramesDir);
        std::filesystem::create_directories(outputVideosDir);
    }

    Menu renderMenu(menuData);

    


    unsigned int iterationCount = 0;
    unsigned int recordedFrameCount = 0;
    double time = 0;
    double dt = DEFAULT_TIME_STEP;
    Solid solid(mesh);

    std::vector<double> timeLog;
    std::vector<double> EnergyLog;
    timeLog.reserve(200);
    EnergyLog.reserve(200);
    timeLog.push_back(time);
    EnergyLog.push_back(0.0);
    while (window.isOpen())
    {
        if(time>14) break;
        PROFILE_SCOPE("Main::frame");
        window.clear(sf::Color::White);
        while (auto event = window.pollEvent())
        {
            ImGui::SFML::ProcessEvent(window, *event);
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }
        
        ImGui::SFML::Update(window, clock.restart());
        ImGui::Text("Data points: %d", (int)timeLog.size());
        renderMenu.processMenu(testRenderer);

        if (renderMenu.menuData.renderInfo.shouldRender)
        {
            
            testRenderer.draw(window);
        }
        if (should_record_structure_evolution)
        {
            const auto windowSize = window.getSize();
            std::vector<std::uint8_t> pixels(windowSize.x * windowSize.y * 4);
            glReadPixels(0, 0,
                         static_cast<GLsizei>(windowSize.x),
                         static_cast<GLsizei>(windowSize.y),
                         GL_RGBA,
                         GL_UNSIGNED_BYTE,
                         pixels.data());

            sf::Image frameImage({windowSize.x, windowSize.y}, pixels.data());
            frameImage.flipVertically();

            const auto framePath = recordingFramesDir / ("frame_" + std::to_string(recordedFrameCount) + ".png");
            if (!frameImage.saveToFile(framePath))
                std::cerr << "Failed to save recording frame: " << framePath << std::endl;
            else
                recordedFrameCount++;
        }
        for(unsigned int i=0; i<5; i++)
        {
            solid.timeStep(dt);
            time += dt;
            timeLog.push_back(time);
            EnergyLog.push_back(solid.computeElasticEnergyFromExtensionCache());
            iterationCount++;
        }
        std::cout << timeLog.size() << " " << EnergyLog.size() << std::endl;
        ImGui::Begin("Energy vs Time");
        ImPlot::SetNextAxesToFit();
        if (ImPlot::BeginPlot("Energy over Time"))
        {
            std::cout<<EnergyLog.data()<<std::endl;
            ImPlot::SetupAxes("time", "energy");
            ImPlot::PlotLine("Energy", timeLog.data(), EnergyLog.data(), static_cast<int>(timeLog.size()));
            ImPlot::EndPlot();
        }
        ImGui::End();
        // std::cout <<"test particle pos "<<mesh.nodePosFlat[0] <<std::endl;
        mesh.syncVecWithFlatVecData();
        testRenderer.updateNodePositionsFromMesh(mesh);
        // profiling::printLiveReport();
        ImGui::SFML::Render(window);
        window.display();
    }
    const std::filesystem::path outputImagesDir{OUTPUT_IMAGES_DIR};
    if (should_record_structure_evolution || shoould_snapshot_energy_log)
    {
        std::filesystem::create_directories(outputImagesDir);
    }

    if (shoould_snapshot_energy_log)
    {
        const auto energyLogCsvPath = outputImagesDir / "energy_log.csv";
        {
            std::ofstream energyLogCsvFile(energyLogCsvPath);
            if (!energyLogCsvFile)
            {
                std::cerr << "Failed to open energy log CSV for writing: " << energyLogCsvPath << std::endl;
            }
            else
            {
                energyLogCsvFile << "time,energy\n";
                for (std::size_t i = 0; i < timeLog.size() && i < EnergyLog.size(); ++i)
                    energyLogCsvFile << timeLog[i] << ',' << EnergyLog[i] << '\n';
            }
        }

        const auto snapshotPath = outputImagesDir / "energy_log_snapshot.png";
        const std::filesystem::path energyPlotScriptPath{SUBROUTINES_PATH};
        const std::string pythonCommand =
            "/bin/python3 \"" + (energyPlotScriptPath / "plot_energy_log.py").string() +
            "\" \"" + energyLogCsvPath.string() + "\" \"" + snapshotPath.string() + "\"";
        std::system(pythonCommand.c_str());
    }

    ImPlot::DestroyContext();

    if (should_record_structure_evolution)
    {
        const auto outputVideoPath = outputVideosDir / "structure_evolution.mp4";
        const std::string ffmpegCommand =
            "ffmpeg -y -framerate 20 -i \"" + recordingFramesDir.string() + "/frame_%d.png\" "
            "-c:v libx264 -pix_fmt yuv420p \"" + outputVideoPath.string() + "\"";
        std::system(ffmpegCommand.c_str());
    }

    return 0;
}