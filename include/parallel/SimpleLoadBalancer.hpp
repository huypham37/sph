#pragma once

#include "parallel/LoadBalancer.hpp"

namespace sph
{
	namespace parallel
	{

		/**
		 * @brief Simple load balancer based on particle count and computation time
		 *
		 * This load balancer adjusts subdomain boundaries to equalize the computational
		 * load across all subdomains. It uses both particle count and the recorded
		 * computation time to make balancing decisions.
		 */
		class SimpleLoadBalancer : public LoadBalancer
		{
		public:
			/**
			 * @brief Constructor
			 *
			 * @param rebalanceInterval Minimum number of steps between rebalancing operations
			 */
			explicit SimpleLoadBalancer(int rebalanceInterval = 30);
			~SimpleLoadBalancer() override = default;

			/**
			 * @brief Determine if rebalancing is needed based on computation times
			 *
			 * @param subdomains Current subdomains with performance metrics
			 * @return true if load imbalance exceeds threshold
			 */
			bool isRebalancingNeeded(const std::vector<std::unique_ptr<Subdomain>> &subdomains) override;

			/**
			 * @brief Rebalance load by adjusting subdomain boundaries
			 *
			 * @param subdomains Vector of all subdomains
			 * @return true if rebalancing was performed
			 */
			bool rebalance(std::vector<std::unique_ptr<Subdomain>> &subdomains) override;

			/**
			 * @brief Set rebalancing interval
			 *
			 * @param steps Minimum number of steps between rebalance operations
			 */
			void setRebalanceInterval(int steps) { rebalanceInterval = steps; }

		private:
			int rebalanceInterval;			 // Minimum steps between rebalancing
			int stepsSinceLastRebalance = 0; // Steps since last rebalance

			/**
			 * @brief Calculate optimal domain boundaries based on load
			 *
			 * @param subdomains Current subdomains
			 * @return Adjusted subdomain boundaries
			 */
			std::vector<std::pair<float, float>> calculateOptimalBoundaries(
				const std::vector<std::unique_ptr<Subdomain>> &subdomains) const;
		};

	} // namespace parallel
} // namespace sph