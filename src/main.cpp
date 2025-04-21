#include <SFML/Graphics.hpp>
#include <iostream>
#include <optional>
#include <string>
#include <sstream>
#include "SPHSimulation.hpp" // Updated to use new class

// Window dimensions
constexpr int WINDOW_WIDTH = 800;
constexpr int WINDOW_HEIGHT = 600;

// Simulation parameters
constexpr float SIMULATION_TIMESTEP = 0.005f;
constexpr int DEFAULT_PARTICLE_COUNT = 2000;
constexpr float MOUSE_FORCE_STRENGTH = 20.0f;

// Particle adjustment parameters
constexpr int PARTICLES_INCREMENT = 100;

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
	simulation.setGravity(0.0f, 9.8f);

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

	// Thread counter text
	sf::Text threadCountText(font, "");
	threadCountText.setCharacterSize(16);
	threadCountText.setFillColor(sf::Color::White);
	threadCountText.setPosition({10, 30});

	// FPS counter text
	sf::Text fpsText(font, "");
	fpsText.setCharacterSize(16);
	fpsText.setFillColor(sf::Color::White);
	fpsText.setPosition({10, 50});

	// Instruction text
	sf::Text instructionsText(font,
							  "Controls:\n"
							  "Move Mouse: Interact with fluid\n"
							  "+/-: Add/Remove 100 particles\n"
							  "[/]: Increase/Decrease threads\n"
							  "P: Toggle parallelization\n"
							  "V: Toggle subdomain visualization\n"
							  "L: Toggle load balancing\n"
							  "B: Toggle load balance visualization\n"
							  "T: Adjust load balance threshold\n"
							  "I: Change load balance interval\n");
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
						simulation.setGravity(0.0f, 9.8f);
					}
					else
					{
						simulation.setGravity(0.0f, -9.8f);
					}
				}

				// Add particles with plus key
				if (keyEvent->code == sf::Keyboard::Key::Add ||
					keyEvent->code == sf::Keyboard::Key::Equal) // = and + share the same key
				{
					simulation.addParticles(PARTICLES_INCREMENT);
				}

				// Remove particles with minus key
				if (keyEvent->code == sf::Keyboard::Key::Subtract ||
					keyEvent->code == sf::Keyboard::Key::Hyphen) // - and _ share the same key
				{
					simulation.removeParticles(PARTICLES_INCREMENT);
				}

				// Increase thread count with ] key
				if (keyEvent->code == sf::Keyboard::Key::RBracket)
				{
					simulation.setNumThreads(simulation.getNumThreads() + 1);
				}

				// Decrease thread count with [ key
				if (keyEvent->code == sf::Keyboard::Key::LBracket)
				{
					simulation.setNumThreads(simulation.getNumThreads() - 1);
				}

				// Toggle parallelization with P key
				if (keyEvent->code == sf::Keyboard::Key::P)
				{
					simulation.setParallelizationEnabled(!simulation.isParallelizationEnabled());
					std::cout << "Parallelization " << (simulation.isParallelizationEnabled() ? "enabled" : "disabled") << std::endl;
				}

				// Toggle subdomain visualization with V key
				if (keyEvent->code == sf::Keyboard::Key::V)
				{
					simulation.setVisualizeSubdomains(!simulation.isVisualizeSubdomains());
					std::cout << "Subdomain visualization " << (simulation.isVisualizeSubdomains() ? "enabled" : "disabled") << std::endl;
				}

				// Toggle load balancing with L key
				if (keyEvent->code == sf::Keyboard::Key::L)
				{
					simulation.setLoadBalancingEnabled(!simulation.isLoadBalancingEnabled());
				}

				// Toggle load balancing visualization with B key
				if (keyEvent->code == sf::Keyboard::Key::B)
				{
					simulation.setVisualizeLoadBalance(!simulation.isVisualizeLoadBalance());
				}

				// Adjust load balance threshold with T key (cycle through 0.1, 0.2, 0.3)
				if (keyEvent->code == sf::Keyboard::Key::T)
				{
					// We don't have a getter for the current threshold, so we'll cycle through values
					static float thresholds[] = {0.1f, 0.2f, 0.3f, 0.5f};
					static int currentThresholdIndex = 0;

					currentThresholdIndex = (currentThresholdIndex + 1) % 4;
					float newThreshold = thresholds[currentThresholdIndex];
					simulation.setLoadBalanceThreshold(newThreshold);
					std::cout << "Load balance threshold set to " << newThreshold << std::endl;
				}

				// Adjust load balance interval with I key (cycle through values)
				if (keyEvent->code == sf::Keyboard::Key::I)
				{
					static int intervals[] = {10, 30, 60, 120};
					static int currentIntervalIndex = 0;

					currentIntervalIndex = (currentIntervalIndex + 1) % 4;
					int newInterval = intervals[currentIntervalIndex];
					simulation.setLoadBalanceInterval(newInterval);
					std::cout << "Load balance interval set to " << newInterval << " frames" << std::endl;
				}
			}
		}

		// Calculate delta time
		float deltaTime = clock.restart().asSeconds();
		accumulator += deltaTime;

		// Update simulation with fixed time step
		while (accumulator >= SIMULATION_TIMESTEP)
		{
			simulation.update(SIMULATION_TIMESTEP);
			accumulator -= SIMULATION_TIMESTEP;
		}

		// Update particle count text
		particleCountText.setString("Particles: " + std::to_string(simulation.getParticleCount()));

		// Update thread count text
		std::stringstream threadText;
		threadText << "Threads: " << simulation.getNumThreads() << "/" << simulation.getMaxThreads();
		if (!simulation.isParallelizationEnabled())
		{
			threadText << " (disabled)";
		}
		threadCountText.setString(threadText.str());

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
		window.draw(threadCountText);
		window.draw(fpsText);
		window.draw(instructionsText);

		// Display what we rendered
		window.display();
	}

	return 0;
}