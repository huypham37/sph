#pragma once

#include "parallel/DomainDecomposer.hpp"
#include <unordered_map>
#include <array>

namespace sph
{
	namespace parallel
	{

		/**
		 * @brief Domain decomposition using Z-order space-filling curve
		 *
		 * This decomposer maps 2D particle positions onto a 1D space-filling curve
		 * (Z-order curve) and then partitions this 1D space into equal segments.
		 * This approach preserves spatial locality better than simple grid decomposition.
		 */
		class SpaceFillingCurveDecomposer : public DomainDecomposer
		{
		public:
			/**
			 * @brief Constructor
			 *
			 * @param resolution Grid resolution for mapping (power of 2, e.g., 1024)
			 */
			explicit SpaceFillingCurveDecomposer(int resolution = 1024);
			~SpaceFillingCurveDecomposer() override = default;

			/**
			 * @brief Create initial domain decomposition using Z-order curve
			 *
			 * @param width Total width of the simulation domain
			 * @param height Total height of the simulation domain
			 * @param numSubdomains Number of subdomains to create
			 * @return Vector of created subdomains
			 */
			std::vector<std::unique_ptr<Subdomain>> createDecomposition(
				float width, float height, int numSubdomains) override;

			/**
			 * @brief Assign particles to subdomains based on Z-order mapping
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
			 * @param subdomains Current subdomains
			 * @param particles All particles in the simulation
			 * @return true if decomposition was changed
			 */
			bool updateDecomposition(
				std::vector<std::unique_ptr<Subdomain>> &subdomains,
				const std::vector<Particle *> &particles) override;

			/**
			 * @brief Set grid resolution for Z-order mapping
			 *
			 * @param res Resolution (must be power of 2)
			 */
			void setResolution(int res);

		private:
			int resolution;		// Grid resolution for Z-order mapping
			float domainWidth;	// Total simulation domain width
			float domainHeight; // Total simulation domain height
			bool needsRebuild;	// Flag indicating if decomposition needs rebuild

			// Cache of subdomain membership for particles
			std::unordered_map<Particle *, int> particleToSubdomainMap;

			// Array of Z-values for each subdomain boundary
			std::vector<uint64_t> subdomainBoundaries;

			/**
			 * @brief Convert 2D position to Z-order index
			 *
			 * @param x X coordinate (0-1 normalized)
			 * @param y Y coordinate (0-1 normalized)
			 * @return Z-order index
			 */
			uint64_t positionToZOrder(float x, float y) const;

			/**
			 * @brief Interleave bits of x and y to create Z-order value
			 *
			 * @param x X coordinate (integer)
			 * @param y Y coordinate (integer)
			 * @return Interleaved Z-order value
			 */
			uint64_t interleave(uint32_t x, uint32_t y) const;

			/**
			 * @brief Find which subdomain a Z-order value belongs to
			 *
			 * @param zValue Z-order value
			 * @return Subdomain index
			 */
			int findSubdomainForZValue(uint64_t zValue) const;
		};

	} // namespace parallel
} // namespace sph