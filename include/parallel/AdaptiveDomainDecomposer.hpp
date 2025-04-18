#pragma once

#include "parallel/DomainDecomposer.hpp"
#include <chrono>

namespace sph
{
	namespace parallel
	{

		/**
		 * @brief Adaptive domain decomposition based on particle density
		 *
		 * This decomposer creates subdomains that adapt to the distribution of particles,
		 * allocating more subdomains in areas with higher particle density to balance load.
		 * It implements a form of Orthogonal Recursive Bisection (ORB) algorithm.
		 */
		class AdaptiveDomainDecomposer : public DomainDecomposer
		{
		public:
			/**
			 * @brief Constructor
			 *
			 * @param adaptivityThreshold Particle density difference required to trigger adaptivity
			 * @param minUpdateInterval Minimum time between updates in milliseconds
			 */
			explicit AdaptiveDomainDecomposer(
				float adaptivityThreshold = 0.2f,
				int minUpdateInterval = 500);

			~AdaptiveDomainDecomposer() override = default;

			/**
			 * @brief Create initial adaptive domain decomposition
			 *
			 * @param width Total width of the simulation domain
			 * @param height Total height of the simulation domain
			 * @param numSubdomains Number of subdomains to create
			 * @return Vector of created subdomains
			 */
			std::vector<std::unique_ptr<Subdomain>> createDecomposition(
				float width, float height, int numSubdomains) override;

			/**
			 * @brief Assign particles to subdomains based on position
			 *
			 * @param particles All particles in the simulation
			 * @param subdomains Vector of subdomains
			 */
			void assignParticlesToSubdomains(
				const std::vector<Particle *> &particles,
				std::vector<std::unique_ptr<Subdomain>> &subdomains) override;

			/**
			 * @brief Update domain decomposition based on particle distribution
			 *
			 * Analyzes current particle distribution and adjusts subdomain boundaries
			 * to better balance computation across subdomains.
			 *
			 * @param subdomains Current subdomains
			 * @param particles All particles in the simulation
			 * @return true if decomposition was changed
			 */
			bool updateDecomposition(
				std::vector<std::unique_ptr<Subdomain>> &subdomains,
				const std::vector<Particle *> &particles) override;

			/**
			 * @brief Set adaptivity threshold
			 *
			 * @param threshold Threshold for when to adapt (0-1)
			 */
			void setAdaptivityThreshold(float threshold) { adaptivityThreshold = threshold; }

			/**
			 * @brief Get current adaptivity threshold
			 *
			 * @return Current threshold value
			 */
			float getAdaptivityThreshold() const { return adaptivityThreshold; }

		private:
			float adaptivityThreshold;							  // Threshold for adapting (between 0-1)
			int minUpdateInterval;								  // Minimum milliseconds between updates
			std::chrono::steady_clock::time_point lastUpdateTime; // Last update time

			/**
			 * @brief Recursively bisect a domain to create subdomains
			 *
			 * @param particles Particles in the current subdomain
			 * @param x Left coordinate of current domain
			 * @param y Top coordinate of current domain
			 * @param width Width of current domain
			 * @param height Height of current domain
			 * @param depth Current recursion depth
			 * @param maxDepth Maximum recursion depth based on numSubdomains
			 * @param subdomains Output vector for created subdomains
			 * @param id Current subdomain ID counter
			 */
			void recursiveBisection(
				const std::vector<Particle *> &particles,
				float x, float y, float width, float height,
				int depth, int maxDepth,
				std::vector<std::unique_ptr<Subdomain>> &subdomains,
				int &id);

			/**
			 * @brief Find best split position for a domain
			 *
			 * @param particles Particles in the domain
			 * @param x Left coordinate
			 * @param y Top coordinate
			 * @param width Width of domain
			 * @param height Height of domain
			 * @param splitVertical Whether to split vertically or horizontally
			 * @return Position of the split
			 */
			float findBestSplitPosition(
				const std::vector<Particle *> &particles,
				float x, float y, float width, float height,
				bool splitVertical);
		};

	} // namespace parallel
} // namespace sph