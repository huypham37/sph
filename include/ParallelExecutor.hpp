#pragma once

#include <vector>
#include <memory>
#include <functional>
#include "Particle.hpp"
#include "parallel/DomainDecomposer.hpp"
#include "parallel/BoundaryManager.hpp"

namespace sph
{

	/**
	 * @brief Manages parallel execution of tasks
	 *
	 * Responsible for delegating work to multiple threads using
	 * domain decomposition, handling ghost regions for boundary particles.
	 */
	class ParallelExecutor
	{
	public:
		/**
		 * @brief Constructor
		 *
		 * @param width Width of simulation domain
		 * @param height Height of simulation domain
		 * @param smoothingRadius SPH smoothing radius for ghost region size
		 */
		ParallelExecutor(float width, float height, float smoothingRadius);

		/**
		 * @brief Destructor
		 */
		~ParallelExecutor();

		/**
		 * @brief Set thread count for parallel execution
		 *
		 * @param count Number of threads to use
		 */
		void setThreadCount(int count);

		/**
		 * @brief Get current thread count
		 *
		 * @return Current number of threads
		 */
		int getThreadCount() const;

		/**
		 * @brief Get maximum available hardware threads
		 *
		 * @return Maximum thread count
		 */
		int getMaxThreadCount() const;

		/**
		 * @brief Set whether parallelization is enabled
		 *
		 * @param enabled True if parallel execution should be used
		 */
		void setParallelizationEnabled(bool enabled);

		/**
		 * @brief Check if parallelization is enabled
		 *
		 * @return True if parallel execution is enabled
		 */
		bool isParallelizationEnabled() const { return parallelizationEnabled; }

		/**
		 * @brief Execute a task in parallel across subdomains
		 *
		 * @param task Function to execute on each subdomain
		 * @param particles Vector of all particles
		 */
		void executeParallel(
			const std::function<void(const std::vector<Particle *> &, const std::vector<Particle *> &)> &task,
			const std::vector<Particle *> &particles);

		/**
		 * @brief Update domain decomposition based on particle positions
		 *
		 * @param particles Vector of all particles
		 */
		void updateDecomposition(const std::vector<Particle *> &particles);

		/**
		 * @brief Get all subdomains for visualization
		 *
		 * @return Reference to subdomains
		 */
		const std::vector<std::unique_ptr<parallel::Subdomain>> &getSubdomains() const { return subdomains; }

		/**
		 * @brief Check if load balancing is needed and perform if necessary
		 *
		 * @return True if load balancing was performed
		 */
		// bool checkAndRebalance();

		/**
		 * @brief Set whether load balancing is enabled
		 *
		 * @param enabled True if load balancing should be used
		 */
		// void setLoadBalancingEnabled(bool enabled) { loadBalancingEnabled = enabled; }

		/**
		 * @brief Check if load balancing is enabled
		 *
		 * @return True if load balancing is enabled
		 */
		// bool isLoadBalancingEnabled() const { return loadBalancingEnabled; }

	private:
		// Domain dimensions
		float width;
		float height;
		float smoothingRadius;

		// Number of threads
		int numThreads;

		// Domain decomposition components
		std::unique_ptr<parallel::DomainDecomposer> domainDecomposer;
		std::unique_ptr<parallel::BoundaryManager> boundaryManager;
		std::vector<std::unique_ptr<parallel::Subdomain>> subdomains;

		bool parallelizationEnabled;

		// Setup procedures
		void initializeComponents();
	};

} // namespace sph