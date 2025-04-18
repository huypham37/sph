#pragma once

#include <vector>
#include <memory>
#include "parallel/Subdomain.hpp"

namespace sph
{
	namespace parallel
	{

		/**
		 * @brief Handles communication and data exchange between subdomains
		 *
		 * Responsible for managing ghost regions and ensuring particles near subdomain
		 * boundaries are properly shared between subdomains.
		 */
		class BoundaryManager
		{
		public:
			/**
			 * @brief Constructor
			 *
			 * @param ghostRegionWidth Width of ghost region around each subdomain (typically SPH smoothing radius)
			 */
			explicit BoundaryManager(float ghostRegionWidth);
			virtual ~BoundaryManager() = default;

			/**
			 * @brief Exchange boundary particles between subdomains
			 *
			 * Identifies particles that need to be shared with neighboring subdomains
			 * and populates ghost particles in each subdomain.
			 *
			 * @param subdomains Vector of all subdomains
			 */
			virtual void exchangeBoundaryData(std::vector<std::unique_ptr<Subdomain>> &subdomains);

			/**
			 * @brief Clear all ghost particles from subdomains
			 *
			 * @param subdomains Vector of all subdomains
			 */
			virtual void clearGhostParticles(std::vector<std::unique_ptr<Subdomain>> &subdomains);

			/**
			 * @brief Set width of ghost regions
			 *
			 * @param width Width of ghost region (should be at least SPH smoothing radius)
			 */
			void setGhostRegionWidth(float width) { ghostRegionWidth = width; }

			/**
			 * @brief Get current width of ghost regions
			 *
			 * @return Ghost region width
			 */
			float getGhostRegionWidth() const { return ghostRegionWidth; }

		protected:
			/**
			 * @brief Determine if two subdomains are neighbors
			 *
			 * @param a First subdomain
			 * @param b Second subdomain
			 * @return true if subdomains share a boundary or corner
			 */
			bool areNeighbors(const Subdomain &a, const Subdomain &b) const;

			/**
			 * @brief Find particles in source subdomain that need to be shared with target
			 *
			 * @param source Source subdomain
			 * @param target Target subdomain
			 * @return Vector of particles that should be ghost particles in target
			 */
			std::vector<Particle *> findParticlesToShare(const Subdomain &source, const Subdomain &target) const;

		private:
			float ghostRegionWidth; // Width of ghost region around each subdomain
		};

	} // namespace parallel
} // namespace sph