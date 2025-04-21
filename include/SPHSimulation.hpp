#pragma once

#include <SFML/Graphics.hpp>
#include <memory>
#include "ParticleSystem.hpp"
#include "SPHPhysics.hpp"
#include "ParallelExecutor.hpp"
#include "Renderer.hpp"

namespace sph
{

	/**
	 * @brief Main coordinator class for the SPH simulation
	 *
	 * Manages the overall simulation, delegating specific responsibilities
	 * to specialized components like physics, particle management,
	 * parallelization, and rendering.
	 */
	class SPHSimulation
	{
	public:
		/**
		 * @brief Constructor
		 *
		 * @param width Width of simulation domain
		 * @param height Height of simulation domain
		 */
		SPHSimulation(float width, float height);

		/**
		 * @brief Destructor
		 */
		~SPHSimulation();

		/**
		 * @brief Update simulation state for one time step
		 *
		 * @param dt Time step
		 */
		void update(float dt);

		/**
		 * @brief Render simulation to SFML window
		 *
		 * @param window SFML window to draw to
		 */
		void draw(sf::RenderWindow &window);

		/**
		 * @brief Initialize simulation with default particle configuration
		 *
		 * @param count Number of particles to create
		 */
		void initializeDefaultParticles(int count);

		/**
		 * @brief Add single particle at specific position
		 *
		 * @param x X-coordinate
		 * @param y Y-coordinate
		 */
		void addParticle(float x, float y);

		/**
		 * @brief Add multiple particles
		 *
		 * @param count Number of particles to add
		 */
		void addParticles(int count);

		/**
		 * @brief Remove specific number of particles
		 *
		 * @param count Number of particles to remove
		 */
		void removeParticles(int count);

		/**
		 * @brief Reset simulation (remove all particles)
		 */
		void reset();

		/**
		 * @brief Apply force to particles near mouse position
		 *
		 * @param mousePos Mouse position in window coordinates
		 * @param strength Force strength (positive = push, negative = pull)
		 */
		void applyMouseForce(const sf::Vector2f &mousePos, float strength);

		// Thread management methods
		void setNumThreads(int num);
		int getNumThreads() const;
		int getMaxThreads() const;

		// Parallelization toggle
		void setParallelizationEnabled(bool enabled);
		bool isParallelizationEnabled() const;

		// Visualization options
		void setVisualizeSubdomains(bool enabled);
		bool isVisualizeSubdomains() const;

		// Simulation parameter setters
		void setGravity(float x, float y);
		void setViscosity(float v);
		void setGasConstant(float k);
		void setRestDensity(float d);

		// Statistics
		size_t getParticleCount() const;

		// Load balancing methods
		void setLoadBalancingEnabled(bool enabled);
		bool isLoadBalancingEnabled() const;
		void setLoadBalanceThreshold(float threshold);
		void setLoadBalanceInterval(int frames);

		// Load balance visualization
		void setVisualizeLoadBalance(bool enabled);
		bool isVisualizeLoadBalance() const;
        
        // Debug options
        void setDebugThreads(bool enabled);
        bool isDebugThreads() const;

	private:
		// Main simulation components
		std::unique_ptr<ParticleSystem> particles;
		std::unique_ptr<SPHPhysics> physics;
		std::unique_ptr<ParallelExecutor> parallelExecutor;
		std::unique_ptr<Renderer> renderer;

		// Simulation dimensions
		float width;
		float height;

		// SPH parameters
		float smoothingRadius;
	};

} // namespace sph