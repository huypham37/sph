#include "parallel/SimpleLoadBalancer.hpp"
#include <algorithm>
#include <numeric>
#include <iostream>

namespace sph
{
	namespace parallel
	{
		SimpleLoadBalancer::SimpleLoadBalancer(int rebalanceInterval)
			: rebalanceInterval(rebalanceInterval)
		{
		}

		bool SimpleLoadBalancer::isRebalancingNeeded(const std::vector<std::unique_ptr<Subdomain>> &subdomains)
		{
			// Don't rebalance if not enough steps have passed since last rebalance
			if (++stepsSinceLastRebalance < rebalanceInterval)
			{
				return false;
			}

			// Need at least 2 subdomains to rebalance
			if (subdomains.size() < 2)
			{
				return false;
			}

			// Calculate average computation time
			float totalTime = 0.0f;
			int totalParticles = 0;

			for (const auto &subdomain : subdomains)
			{
				totalTime += subdomain->getLastComputationTime();
				totalParticles += subdomain->getParticles().size();
			}

			// Avoid division by zero
			if (totalTime <= 0.0f || totalParticles <= 0)
			{
				return false;
			}

			float avgTime = totalTime / subdomains.size();

			// Check for imbalance
			for (const auto &subdomain : subdomains)
			{
				float time = subdomain->getLastComputationTime();
				// If any subdomain is significantly slower than average, rebalance
				if (time > avgTime * (1.0f + imbalanceThreshold))
				{
					stepsSinceLastRebalance = 0;
					std::cout << "Load imbalance detected: " << time << " vs avg " << avgTime << std::endl;
					return true;
				}
			}

			return false;
		}

		bool SimpleLoadBalancer::rebalance(std::vector<std::unique_ptr<Subdomain>> &subdomains)
		{
			if (subdomains.size() < 2)
			{
				return false;
			}

			// For now, implement a simple 1D strip decomposition adjustment
			// This assumes subdomains are arranged horizontally

			// Calculate optimal boundaries
			auto newBoundaries = calculateOptimalBoundaries(subdomains);

			// Apply new boundaries (in a real implementation, this would be more sophisticated)
			for (size_t i = 0; i < subdomains.size(); ++i)
			{
				// In a real implementation, you would redistribute particles based on new boundaries
				// For this simple version, we're just noting that rebalancing occurred
				std::cout << "Rebalanced subdomain " << i << " boundary: "
						  << newBoundaries[i].first << " -> " << newBoundaries[i].second << std::endl;
			}

			// In a real implementation, you would now need to:
			// 1. Reassign particles to their new subdomains
			// 2. Update ghost particles
			// 3. Clear performance metrics to start fresh measurements

			return true;
		}

		std::vector<std::pair<float, float>> SimpleLoadBalancer::calculateOptimalBoundaries(
			const std::vector<std::unique_ptr<Subdomain>> &subdomains) const
		{
			// For this simple implementation, we'll calculate boundaries proportionally to computation time

			std::vector<std::pair<float, float>> boundaries(subdomains.size());

			// Extract computation times
			std::vector<float> times;
			float totalTime = 0.0f;

			for (const auto &subdomain : subdomains)
			{
				float time = subdomain->getLastComputationTime();
				// If time is too small, use particle count as a proxy
				if (time < 0.001f)
				{
					time = static_cast<float>(subdomain->getParticles().size());
				}
				times.push_back(time);
				totalTime += time;
			}

			// Avoid division by zero
			if (totalTime <= 0.0f)
			{
				// Return original boundaries if we can't calculate new ones
				for (size_t i = 0; i < subdomains.size(); ++i)
				{
					boundaries[i] = {subdomains[i]->getX(), subdomains[i]->getX() + subdomains[i]->getWidth()};
				}
				return boundaries;
			}

			// Calculate total domain width
			float minX = subdomains[0]->getX();
			float maxX = subdomains[0]->getX() + subdomains[0]->getWidth();

			for (size_t i = 1; i < subdomains.size(); ++i)
			{
				minX = std::min(minX, subdomains[i]->getX());
				maxX = std::max(maxX, subdomains[i]->getX() + subdomains[i]->getWidth());
			}

			float totalWidth = maxX - minX;

			// Calculate proportional boundaries
			float currentX = minX;
			for (size_t i = 0; i < subdomains.size(); ++i)
			{
				float proportion = times[i] / totalTime;
				float width = totalWidth * proportion;

				boundaries[i] = {currentX, currentX + width};
				currentX += width;
			}

			// Ensure the last boundary ends exactly at maxX to avoid floating-point errors
			boundaries.back().second = maxX;

			return boundaries;
		}

	} // namespace parallel
} // namespace sph