#include <SFML/Graphics.hpp>
#include <iostream>
#include <optional>
#include <string>
#include <sstream>
#include "SPHSimulation.hpp" // Updated to use new class
#include "SPHConfig.hpp"

// Window dimensions
constexpr int WINDOW_WIDTH = Config::WIDTH;
constexpr int WINDOW_HEIGHT = Config::HEIGHT;

// Simulation parameters
constexpr float SIMULATION_TIMESTEP = Config::TIME_STEP;
constexpr int DEFAULT_PARTICLE_COUNT = Config::PARTICLE_COUNT;
constexpr float MOUSE_FORCE_STRENGTH = 20.0f;
constexpr float GRAVITY_Y = Config::GRAVITY;

// FPS counter parameters
constexpr float FPS_UPDATE_INTERVAL = 0.5f; // Update FPS display every 0.5 seconds

int main()
{
	// Create the main window
	sf::RenderWindow window(sf::VideoMode({WINDOW_WIDTH, WINDOW_HEIGHT}), "SPH Fluid Simulation", sf::Style::Close);
	window.setFramerateLimit(60);

	// Create UI font
	sf::Font font;
	if (!font.openFromFile("/System/Library/Fonts/Helvetica.ttc")) // Using openFromFile for SFML 3.0
	{
		std::cerr << "Error loading font" << std::endl;
		return 1;
	}

	// Create SPH simulation
	sph::SPHSimulation simulation(WINDOW_WIDTH, WINDOW_HEIGHT); // Updated to use new class

	// Initialize gravity
	simulation.setGravity(0.0f, GRAVITY_Y);

	// Initialize the simulation with default particles
	simulation.initializeDefaultParticles(DEFAULT_PARTICLE_COUNT);

	// Timer for simulation steps
	sf::Clock clock;
	float accumulator = 0.0f;

	// Particle counter text
	sf::Text particleCountText(font, "");
	particleCountText.setCharacterSize(16);
	particleCountText.setFillColor(sf::Color::White);
	particleCountText.setPosition({10, 10});

	// FPS counter text
	sf::Text fpsText(font, "");
	fpsText.setCharacterSize(16);
	fpsText.setFillColor(sf::Color::White);
	fpsText.setPosition({10, 30});

	// Instruction text
	sf::Text instructionsText(font,
							  "Controls:\n"
							  "Space: Pause/Resume simulation\n"
							  "Move Mouse: Interact with fluid\n");
	instructionsText.setCharacterSize(14);
	instructionsText.setFillColor(sf::Color::White);
	instructionsText.setPosition({10, WINDOW_HEIGHT - 150});

	// FPS counter variables
	sf::Clock fpsClock;
	int frameCount = 0;

	// Main game loop
	bool gravityDown = true;
	sf::Vector2f lastMousePos;
	bool isMouseInWindow = false;
	bool simulationPaused = false;

	while (window.isOpen())
	{
		// Process events
		std::optional<sf::Event> eventOpt;
		while ((eventOpt = window.pollEvent()))
		{
			const auto &event = *eventOpt;

			if (event.is<sf::Event::Closed>())
			{
				window.close();
			}

			// Mouse enter/leave window tracking
			if (event.is<sf::Event::MouseEntered>())
			{
				isMouseInWindow = true;
			}
			else if (event.is<sf::Event::MouseLeft>())
			{
				isMouseInWindow = false;
			}

			// Mouse movement for fluid interaction
			if (event.is<sf::Event::MouseMoved>())
			{
				if (isMouseInWindow)
				{
					const auto *moveEvent = event.getIf<sf::Event::MouseMoved>();
					// Fixed: Access position vector instead of x/y members
					sf::Vector2f currentMousePos = {
						static_cast<float>(moveEvent->position.x),
						static_cast<float>(moveEvent->position.y)};

					// Apply force to particles based on mouse movement
					simulation.applyMouseForce(currentMousePos, MOUSE_FORCE_STRENGTH);

					// Save current position for next frame
					lastMousePos = currentMousePos;
				}
			}

			// Keyboard input
			if (event.is<sf::Event::KeyPressed>())
			{
				const auto *keyEvent = event.getIf<sf::Event::KeyPressed>();

				// Reset simulation with R key
				if (keyEvent->code == sf::Keyboard::Key::R)
				{
					simulation.initializeDefaultParticles(DEFAULT_PARTICLE_COUNT);
				}

				// Toggle gravity direction with G key
				if (keyEvent->code == sf::Keyboard::Key::G)
				{
					gravityDown = !gravityDown;
					if (gravityDown)
					{
						simulation.setGravity(0.0f, GRAVITY_Y);
					}
					else
					{
						simulation.setGravity(0.0f, -GRAVITY_Y);
					}
				}

				// Toggle pause with space bar
				if (keyEvent->code == sf::Keyboard::Key::Space)
				{
					simulationPaused = !simulationPaused;
					std::cout << "Simulation " << (simulationPaused ? "paused" : "resumed") << std::endl;
				}
			}
		}

		// Calculate delta time
		float deltaTime = clock.restart().asSeconds();
		accumulator += deltaTime;

		// Update simulation with fixed time step
		while (accumulator >= SIMULATION_TIMESTEP && !simulationPaused)
		{
			simulation.update(SIMULATION_TIMESTEP);
			accumulator -= SIMULATION_TIMESTEP;
		}

		// Update particle count text
		particleCountText.setString("Particles: " + std::to_string(simulation.getParticleCount()));

		// Update FPS counter
		frameCount++;
		if (fpsClock.getElapsedTime().asSeconds() >= FPS_UPDATE_INTERVAL)
		{
			float fps = frameCount / fpsClock.getElapsedTime().asSeconds();
			fpsText.setString("FPS: " + std::to_string(static_cast<int>(fps)));
			fpsClock.restart();
			frameCount = 0;
		}

		// Clear the window
		window.clear(sf::Color(30, 30, 40));

		// Draw simulation
		simulation.draw(window);

		// Draw UI
		window.draw(particleCountText);
		window.draw(fpsText);
		window.draw(instructionsText);

		// Display what we rendered
		window.display();
	}

	return 0;
}