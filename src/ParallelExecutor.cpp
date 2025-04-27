#include "ParallelExecutor.hpp"
#include <omp.h>
#include <iostream>
#include <algorithm>
#include "parallel/GridDomainDecomposer.hpp"
#include "parallel/SimpleLoadBalancer.hpp"

namespace sph
{

	ParallelExecutor::ParallelExecutor(float width, float height, float smoothingRadius)
		: width(width),
		  height(height),
		  smoothingRadius(smoothingRadius),
		  parallelizationEnabled(true),
		  loadBalancingEnabled(true),
		  numThreads(omp_get_max_threads())
	{
		// Initialize OpenMP with maximum available threads
		omp_set_num_threads(numThreads);
		std::cout << "Parallel executor initialized with " << numThreads << " threads" << std::endl;

		// Create components for domain decomposition
		initializeComponents();
	}

	ParallelExecutor::~ParallelExecutor()
	{
		// Unique pointers will be automatically cleaned up
	}

	void ParallelExecutor::initializeComponents()
	{
		// Create GridDomainDecomposer by default
		domainDecomposer = std::make_unique<parallel::GridDomainDecomposer>();

		// Create BoundaryManager with smoothing radius as ghost region width
		boundaryManager = std::make_unique<parallel::BoundaryManager>(smoothingRadius);

		// Create SimpleLoadBalancer with default settings
		loadBalancer = std::make_unique<parallel::SimpleLoadBalancer>();
		std::cout << "Load balancer initialized with default settings" << std::endl;
	}

	void ParallelExecutor::setThreadCount(int count)
	{
		// Ensure thread count is at least 1 and at most the hardware maximum
		int maxThreads = getMaxThreadCount();
		count = std::max(1, std::min(count, maxThreads));

		if (count != numThreads)
		{
			numThreads = count;

			// Update OpenMP thread count
			omp_set_num_threads(numThreads);
			std::cout << "Thread count set to " << numThreads << " of " << maxThreads << " available" << std::endl;

			// Force recreation of subdomains with new thread count
			if (parallelizationEnabled)
			{
				subdomains.clear();
			}
		}
	}

	int ParallelExecutor::getThreadCount() const
	{
		return numThreads;
	}

	int ParallelExecutor::getMaxThreadCount() const
	{
		return omp_get_max_threads();
	}

	void ParallelExecutor::setParallelizationEnabled(bool enabled)
	{
		if (parallelizationEnabled != enabled)
		{
			parallelizationEnabled = enabled;
			std::cout << "Parallelization " << (enabled ? "enabled" : "disabled") << std::endl;
		}
	}

	void ParallelExecutor::updateDecomposition(const std::vector<Particle *> &particles)
	{
		if (!parallelizationEnabled || particles.empty())
		{
			return;
		}

		// Create subdomains if they don't exist yet
		if (subdomains.empty())
		{
			subdomains = domainDecomposer->createDecomposition(width, height, numThreads * 3);
		}

		// Assign particles to subdomains
		domainDecomposer->assignParticlesToSubdomains(particles, subdomains);

		// Exchange boundary data between subdomains
		boundaryManager->exchangeBoundaryData(subdomains);
	}

	bool ParallelExecutor::checkAndRebalance()
	{
		if (!parallelizationEnabled || !loadBalancingEnabled || subdomains.empty())
		{
			return false;
		}

		// Check if rebalancing is needed
		if (loadBalancer->isRebalancingNeeded(subdomains))
		{
			std::cout << "Load balancing in progress..." << std::endl;
			bool rebalanced = loadBalancer->rebalance(subdomains);
			if (rebalanced)
			{
				std::cout << "Load balancing completed" << std::endl;
			}
			return rebalanced;
		}
		return false;
	}

	void ParallelExecutor::setLoadBalanceThreshold(float threshold)
	{
		if (loadBalancer)
		{
			loadBalancer->setImbalanceThreshold(threshold);
			std::cout << "Load balance threshold set to " << threshold << std::endl;
		}
	}

	void ParallelExecutor::setLoadBalanceInterval(int frames)
	{
		if (auto simpleBalancer = dynamic_cast<parallel::SimpleLoadBalancer *>(loadBalancer.get()))
		{
			simpleBalancer->setRebalanceInterval(frames);
			std::cout << "Load balance interval set to " << frames << " frames" << std::endl;
		}
	}

	void ParallelExecutor::executeParallel(
		const std::function<void(const std::vector<Particle *> &, const std::vector<Particle *> &)> &task,
		const std::vector<Particle *> &particles)
	{
		if (!parallelizationEnabled || subdomains.empty())
		{
			// If parallelization is disabled, execute the task sequentially with all particles
			task(particles, {});
			return;
		}

		// Process each subdomain in parallel using OpenMP
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < static_cast<int>(subdomains.size()); i++)
		{
			auto &subdomain = subdomains[i];

			const auto &subdomainParticles = subdomain->getParticles();
			const auto &ghostParticles = subdomain->getGhostParticles();

			// Execute the task on this subdomain's particles
			task(subdomainParticles, ghostParticles);
		}

		// Check if load balancing is needed after computation
		if (loadBalancingEnabled)
		{
			checkAndRebalance();
		}
	}

} // namespace sph