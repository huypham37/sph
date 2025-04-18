#pragma once

#include <vector>
#include <memory>
#include "parallel/Subdomain.hpp"

namespace sph
{
	namespace parallel
	{

		/**
		 * @brief Interface for load balancing strategies
		 *
		 * This abstract class defines the interface for different load balancing approaches
		 * to distribute computational workloads evenly across subdomains.
		 */
		class LoadBalancer
		{
		public:
			virtual ~LoadBalancer() = default;

			/**
			 * @brief Analyze current load distribution and determine if rebalancing is needed
			 *
			 * @param subdomains Current subdomains with performance metrics
			 * @return true if load imbalance exceeds threshold and rebalancing is needed
			 */
			virtual bool isRebalancingNeeded(const std::vector<std::unique_ptr<Subdomain>> &subdomains) = 0;

			/**
			 * @brief Rebalance workload across subdomains
			 *
			 * @param subdomains Vector of all subdomains
			 * @return true if rebalancing was performed, false otherwise
			 */
			virtual bool rebalance(std::vector<std::unique_ptr<Subdomain>> &subdomains) = 0;

			/**
			 * @brief Set the imbalance threshold
			 *
			 * Sets the threshold at which load imbalance triggers rebalancing.
			 * For example, 0.2 means rebalance if any subdomain's time is more than 20%
			 * above the average.
			 *
			 * @param threshold Value between 0 and 1
			 */
			void setImbalanceThreshold(float threshold) { imbalanceThreshold = threshold; }

			/**
			 * @brief Get current imbalance threshold
			 *
			 * @return Current threshold value
			 */
			float getImbalanceThreshold() const { return imbalanceThreshold; }

		protected:
			float imbalanceThreshold = 0.2f; // Default 20% threshold
		};

	} // namespace parallel
} // namespace sph