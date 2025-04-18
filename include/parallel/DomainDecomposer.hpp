#pragma once

#include <vector>
#include <memory>
#include "parallel/Subdomain.hpp"
#include "Particle.hpp"

namespace sph
{
	namespace parallel
	{

		/**
		 * @brief Interface for domain decomposition strategies
		 *
		 * This abstract class defines the interface for different domain
		 * decomposition approaches (grid-based, octree, space-filling curves, etc.)
		 */
		class DomainDecomposer
		{
		public:
			virtual ~DomainDecomposer() = default;

			/**
			 * @brief Create initial domain decomposition
			 *
			 * @param width Total width of the simulation domain
			 * @param height Total height of the simulation domain
			 * @param numSubdomains Number of subdomains to create
			 * @return Vector of created subdomains
			 */
			virtual std::vector<std::unique_ptr<Subdomain>> createDecomposition(
				float width, float height, int numSubdomains) = 0;

			/**
			 * @brief Assign particles to appropriate subdomains
			 *
			 * @param particles All particles in the simulation
			 * @param subdomains Vector of subdomains
			 */
			virtual void assignParticlesToSubdomains(
				const std::vector<Particle *> &particles,
				std::vector<std::unique_ptr<Subdomain>> &subdomains) = 0;

			/**
			 * @brief Update decomposition based on current particle distribution
			 *
			 * Used for dynamic/adaptive domain decomposition strategies.
			 *
			 * @param subdomains Current subdomains
			 * @param particles All particles in the simulation
			 * @return true if decomposition was changed, false otherwise
			 */
			virtual bool updateDecomposition(
				std::vector<std::unique_ptr<Subdomain>> &subdomains,
				const std::vector<Particle *> &particles) = 0;
		};

	} // namespace parallel
} // namespace sph