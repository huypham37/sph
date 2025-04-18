#pragma once

#include "parallel/DomainDecomposer.hpp"

namespace sph
{
	namespace parallel
	{

		/**
		 * @brief Simple grid-based domain decomposition strategy
		 *
		 * Decomposes the simulation domain into a regular grid of subdomains.
		 * This is the simplest domain decomposition approach with minimal overhead.
		 */
		class GridDomainDecomposer : public DomainDecomposer
		{
		public:
			/**
			 * @brief Constructor
			 */
			GridDomainDecomposer() = default;
			~GridDomainDecomposer() override = default;

			/**
			 * @brief Create initial grid-based domain decomposition
			 *
			 * @param width Total width of the simulation domain
			 * @param height Total height of the simulation domain
			 * @param numSubdomains Number of subdomains to create
			 * @return Vector of created subdomains
			 */
			std::vector<std::unique_ptr<Subdomain>> createDecomposition(
				float width, float height, int numSubdomains) override;

			/**
			 * @brief Assign particles to appropriate grid cells
			 *
			 * @param particles All particles in the simulation
			 * @param subdomains Vector of subdomains
			 */
			void assignParticlesToSubdomains(
				const std::vector<Particle *> &particles,
				std::vector<std::unique_ptr<Subdomain>> &subdomains) override;

			/**
			 * @brief Update decomposition based on current particle distribution
			 *
			 * In static grid decomposition, this is a no-op since grid boundaries don't change
			 *
			 * @param subdomains Current subdomains
			 * @param particles All particles in the simulation
			 * @return Always returns false since grid doesn't change
			 */
			bool updateDecomposition(
				std::vector<std::unique_ptr<Subdomain>> &subdomains,
				const std::vector<Particle *> &particles) override;

		private:
			/**
			 * @brief Calculate optimal grid dimensions for given number of subdomains
			 *
			 * @param numSubdomains Target number of subdomains
			 * @param width Domain width
			 * @param height Domain height
			 * @param gridCols Output number of columns
			 * @param gridRows Output number of rows
			 */
			void calculateGridDimensions(
				int numSubdomains, float width, float height,
				int &gridCols, int &gridRows) const;
		};

	} // namespace parallel
} // namespace sph