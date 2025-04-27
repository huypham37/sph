#pragma once

#include <vector>
#include "Particle.hpp"
#include "Grid.hpp"

namespace sph
{

	/**
	 * @brief Manages system of particles and spatial partitioning
	 *
	 * Responsible for managing all particles in the simulation,
	 * their creation, deletion, and spatial partitioning with Grid.
	 */
	class ParticleSystem
	{
	public:
		/**
		 * @brief Constructor
		 *
		 * @param width Width of simulation domain
		 * @param height Height of simulation domain
		 * @param smoothingRadius SPH smoothing radius
		 */
		ParticleSystem(float width, float height, float smoothingRadius);

		/**
		 * @brief Destructor
		 */
		~ParticleSystem();

		/**
		 * @brief Add particle at specific position
		 *
		 * @param x X-coordinate
		 * @param y Y-coordinate
		 * @return Pointer to the created particle
		 */
		Particle *addParticle(float x, float y);

		/**
		 * @brief Add multiple particles in a small area
		 *
		 * @param count Number of particles to add
		 */
		void addParticles(int count);

		/**
		 * @brief Remove a specific number of particles
		 *
		 * @param count Number of particles to remove
		 */
		void removeParticles(int count);

		/**
		 * @brief Initialize default particle configuration
		 *
		 * @param count Number of particles to create
		 */
		void initialize(int count);

		/**
		 * @brief Initialize particles in a dam break configuration
		 *
		 * Creates a dense block of particles in one side of the domain
		 * that will collapse and flow when simulation starts
		 *
		 * @param count Number of particles to create
		 * @param damWidth Width of the dam as a fraction of domain width (0.0-1.0)
		 * @param damHeight Height of the dam as a fraction of domain height (0.0-1.0)
		 */
		void initializeDamBreak(int count, float damWidth = 0.4f, float damHeight = 0.8f);

		/**
		 * @brief Delete all particles
		 */
		void reset();

		/**
		 * @brief Update grid for spatial partitioning
		 */
		void updateGrid();

		/**
		 * @brief Get reference to all particles
		 *
		 * @return Reference to particles vector
		 */
		const std::vector<Particle *> &getParticles() const { return particles; }

		/**
		 * @brief Get pointer to spatial partitioning grid
		 *
		 * @return Pointer to Grid
		 */
		Grid *getGrid() const { return grid; }

		/**
		 * @brief Get total number of particles
		 *
		 * @return Particle count
		 */
		size_t getParticleCount() const { return particles.size(); }

		/**
		 * @brief Get simulation domain dimensions
		 *
		 * @return Width of simulation domain
		 */
		float getWidth() const { return width; }

		/**
		 * @brief Get simulation domain dimensions
		 *
		 * @return Height of simulation domain
		 */
		float getHeight() const { return height; }

	private:
		std::vector<Particle *> particles; // All particles in simulation
		Grid *grid;						   // Spatial partitioning grid
		float width;					   // Width of simulation domain
		float height;					   // Height of simulation domain
		float smoothingRadius;			   // SPH smoothing radius
	};

} // namespace sph