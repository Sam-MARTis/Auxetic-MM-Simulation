#include <iostream>
#include <SFML/Graphics.hpp>
#include "renderer.hpp"
#include "utils.hpp"


int main()
{
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
    RendererSetupInformation setupInfo;
    setupInfo.renderOrigin = Vec2D(DEFAULT_SCREEN_WIDTH/2.0, DEFAULT_SCREEN_HEIGHT/2.0);
    setupInfo.nodePositions = particlePositions;
    Renderer testRenderer { setupInfo };


    sf::ContextSettings contextSettings;
    contextSettings.antiAliasingLevel = 16;
    sf::RenderWindow window = sf::RenderWindow(sf::VideoMode({static_cast<unsigned int>(round(DEFAULT_SCREEN_WIDTH)), static_cast<unsigned int>(round(DEFAULT_SCREEN_HEIGHT))}),"Testing");
    window.setFramerateLimit(60);
    while (window.isOpen())
    {
        window.clear(sf::Color::Black);
        while (auto event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }
        testRenderer.draw(window);
        window.display();
    }
    return 0;
}