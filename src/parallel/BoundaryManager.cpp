#include "parallel/BoundaryManager.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>

namespace sph
{
	namespace parallel
	{
		BoundaryManager::BoundaryManager(float ghostRegionWidth)
			: ghostRegionWidth(ghostRegionWidth)
		{
		}

		void BoundaryManager::exchangeBoundaryData(std::vector<std::unique_ptr<Subdomain>> &subdomains)
		{
			// Clear all ghost particles first
			clearGhostParticles(subdomains);

			// Update boundary regions for all subdomains
			for (auto &subdomain : subdomains)
			{
				subdomain->updateBoundaryRegions(ghostRegionWidth);
			}

			// For each pair of subdomains, exchange particles near boundaries
			for (size_t i = 0; i < subdomains.size(); ++i)
			{
				auto &source = subdomains[i];

				for (size_t j = 0; j < subdomains.size(); ++j)
				{
					// Skip self
					if (i == j)
						continue;

					auto &target = subdomains[j];

					// Check if these subdomains are neighbors
					if (areNeighbors(*source, *target))
					{
						// Time measurement for finding particles to share
						auto startTime = std::chrono::high_resolution_clock::now();
						
						// Find particles in source subdomain that need to be shared with target
						auto particlesToShare = findParticlesToShare(*source, *target);
						
						auto endTime = std::chrono::high_resolution_clock::now();
						auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
						

						// Add these particles as ghost particles to the target subdomain
						for (auto *particle : particlesToShare)
						{
							target->addGhostParticle(particle);
						}
					}
				}
			}
		}

		void BoundaryManager::clearGhostParticles(std::vector<std::unique_ptr<Subdomain>> &subdomains)
		{
			for (auto &subdomain : subdomains)
			{
				subdomain->clearGhostParticles();
			}
		}

		bool BoundaryManager::areNeighbors(const Subdomain &a, const Subdomain &b) const
		{
			// Define the extended bounds to account for diagonal neighbors as well
			float ax1 = a.getX() - ghostRegionWidth;
			float ay1 = a.getY() - ghostRegionWidth;
			float ax2 = a.getX() + a.getWidth() + ghostRegionWidth;
			float ay2 = a.getY() + a.getHeight() + ghostRegionWidth;

			float bx1 = b.getX();
			float by1 = b.getY();
			float bx2 = b.getX() + b.getWidth();
			float by2 = b.getY() + b.getHeight();

			// Check if the extended bounds of subdomain A overlap with subdomain B
			bool overlapsX = (bx2 > ax1) && (bx1 < ax2);
			bool overlapsY = (by2 > ay1) && (by1 < ay2);

			return overlapsX && overlapsY;
		}

		std::vector<Particle *> BoundaryManager::findParticlesToShare(const Subdomain &source, const Subdomain &target) const
		{
			std::vector<Particle *> particlesToShare;

			// If subdomains aren't neighbors, return empty vector
			if (!areNeighbors(source, target))
			{
				return particlesToShare;
			}

			// Get shared boundary regions between source and target
			auto sharedRegions = source.getSharedBoundaryRegions(target, ghostRegionWidth);

			// For each shared region, check particles in that region
			for (const auto &region : sharedRegions)
			{
				const auto &regionParticles = source.getParticlesInBoundaryRegion(region);

				// Get extended bounds of target subdomain
				float tx1 = target.getX() - ghostRegionWidth;
				float ty1 = target.getY() - ghostRegionWidth;
				float tx2 = target.getX() + target.getWidth() + ghostRegionWidth;
				float ty2 = target.getY() + target.getHeight() + ghostRegionWidth;

				// Only add particles that are actually within the extended target boundary
				for (const auto *particle : regionParticles)
				{
					sf::Vector2f pos = particle->getPosition();

					// If particle is within ghost region of target, add it
					if (pos.x >= tx1 && pos.x < tx2 && pos.y >= ty1 && pos.y < ty2)
					{
						particlesToShare.push_back(const_cast<Particle *>(particle));
					}
				}
			}

			return particlesToShare;
		}
	}
}