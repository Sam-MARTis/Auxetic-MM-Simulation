#include <iostream>
#include <SFML/Graphics.hpp>

#include "imgui-SFML.h"

#include "utils.hpp"
#include "renderer.hpp"
#include "menu.hpp"

int main()
{
    sf::Clock clock;
    std::cout<<"Hello World!"<<std::endl;
    unsigned int numParticles = 100;
    std::vector<Vec2D> particlePositions;
    particlePositions.reserve(numParticles);
    for (int i = 0; i < numParticles; i++)
    {
        Vec2D particle;
        particle.x = randValueMine(-1.0, 1.0);
        particle.y = randValueMine(-1.0, 1.0);
        particlePositions.push_back(particle);
    }
    RendererSetupInformation setupInfo{.nodePositions = particlePositions};
    setupInfo.renderOrigin = Vec2D(DEFAULT_SCREEN_WIDTH/2.0, DEFAULT_SCREEN_HEIGHT/2.0);
    RendererUpdateData renderUpdateInfo{
        .propertyUpdated = RenderUpdateFlags::NONE,
        .info = setupInfo
    };
    SimulationData simInfo{
    .shouldTimeStep = true,
    .shouldRender = true,
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
    ImGui::SFML::Init(window);
    window.setFramerateLimit(60);

    Menu renderMenu(menuData);

    while (window.isOpen())
    {
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
        testRenderer.draw(window);
        ImGui::SFML::Render(window);
        window.display();
    }
    return 0;
}