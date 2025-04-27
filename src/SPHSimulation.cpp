#include "SPHSimulation.hpp"
#include <iostream>
#include <chrono>

namespace sph
{

	SPHSimulation::SPHSimulation(float width, float height)
		: width(width),
		  height(height),
		  smoothingRadius(14.0f)
	{
		// Create component instances
		particles = std::make_unique<ParticleSystem>(width, height, smoothingRadius);
		physics = std::make_unique<SPHPhysics>();
		physics->setSmoothingRadius(smoothingRadius);
		parallelExecutor = std::make_unique<ParallelExecutor>(width, height, smoothingRadius);
		renderer = std::make_unique<Renderer>();

		std::cout << "SPH Simulation initialized with dimensions " << width << "x" << height << std::endl;
	}

	SPHSimulation::~SPHSimulation()
	{
		// Unique pointers will clean up automatically
	}

	void SPHSimulation::update(float dt)
	{
		constexpr int sub_steps = 1; // Minimal substeps for stability
		float sub_dt = dt / sub_steps;

		// Update grid and cache neighbors once per frame
		particles->updateGrid();
		for (auto *particle : particles->getParticles())
		{
			particle->cachedNeighbors = particles->getGrid()->getNeighbors(particle, smoothingRadius * 2);
		}

		// Substep loop for physics
		for (int i = 0; i < sub_steps; ++i)
		{
			if (parallelExecutor->isParallelizationEnabled() && particles->getParticleCount() > 1000) // Threshold for parallelization
			{
				// Density and pressure task
				auto densityPressureTask = [this](const std::vector<Particle *> &localParticles,
												  const std::vector<Particle *> &)
				{
					physics->computeDensityPressure(localParticles, particles->getGrid());
				};

				auto forcesTask = [this](const std::vector<Particle *> &localParticles,
                                         const std::vector<Particle *> &)
                {
                    physics->computeForces(localParticles, particles->getGrid());
                };

				auto boundaryForcesTask = [this](const std::vector<Particle *> &localParticles,
					const std::vector<Particle *> &)
				{
					physics->computeBoundaryForces(localParticles, width, height);
				};

				// Integration task
				auto integrateTask = [this, sub_dt](const std::vector<Particle *> &localParticles,
													const std::vector<Particle *> &)
				{
					physics->integrate(localParticles, sub_dt);
				};

				// Execute in parallel
				parallelExecutor->executeParallel(densityPressureTask, particles->getParticles());
				parallelExecutor->executeParallel(forcesTask, particles->getParticles());
				parallelExecutor->executeParallel(boundaryForcesTask, particles->getParticles());
				parallelExecutor->executeParallel(integrateTask, particles->getParticles());
			}
			else
			{
				// Sequential execution for small particle counts
				physics->computeDensityPressure(particles->getParticles(), particles->getGrid());
				physics->computeForces(particles->getParticles(), particles->getGrid());
				physics->computeBoundaryForces(particles->getParticles(), width, height);
				physics->integrate(particles->getParticles(), sub_dt);
			}

			// Update particle positions
			for (auto *particle : particles->getParticles())
			{
				particle->update(sub_dt);
			}
		}

		// Resolve collisions once per frame
		physics->resolveCollisions(particles->getParticles(), particles->getGrid(), width, height);
	}

	void SPHSimulation::draw(sf::RenderWindow &window)
	{
		// Only update visuals for particles that are visible
		// You might need a culling step to determine which particles are on-screen
		for (auto *particle : particles->getParticles())
		{
			particle->updateVisuals();
		}

		// Draw particles
		renderer->drawParticles(particles->getParticles(), window);

		// Draw subdomains if enabled
		if (renderer->isVisualizeSubdomains() && parallelExecutor->isParallelizationEnabled())
		{
			renderer->drawSubdomains(parallelExecutor->getSubdomains(), window);
		}

		// Draw load balance info if enabled
		if (renderer->isVisualizeLoadBalance() && parallelExecutor->isParallelizationEnabled())
		{
			renderer->drawSubdomainLoadInfo(parallelExecutor->getSubdomains(), window);
		}
	}

	void SPHSimulation::addParticle(float x, float y)
	{
		particles->addParticle(x, y);
	}

	void SPHSimulation::addParticles(int count)
	{
		particles->addParticles(count);
	}

	void SPHSimulation::removeParticles(int count)
	{
		particles->removeParticles(count);
	}

	void SPHSimulation::reset()
	{
		particles->reset();
	}

	void SPHSimulation::initializeDefaultParticles(int count)
	{
		particles->initializeDamBreak(count);
	}

	void SPHSimulation::setNumThreads(int num)
	{
		parallelExecutor->setThreadCount(num);
	}

	int SPHSimulation::getNumThreads() const
	{
		return parallelExecutor->getThreadCount();
	}

	int SPHSimulation::getMaxThreads() const
	{
		return parallelExecutor->getMaxThreadCount();
	}

	void SPHSimulation::setParallelizationEnabled(bool enabled)
	{
		parallelExecutor->setParallelizationEnabled(enabled);
	}

	bool SPHSimulation::isParallelizationEnabled() const
	{
		return parallelExecutor->isParallelizationEnabled();
	}

	void SPHSimulation::setVisualizeSubdomains(bool enabled)
	{
		renderer->setVisualizeSubdomains(enabled);
	}

	bool SPHSimulation::isVisualizeSubdomains() const
	{
		return renderer->isVisualizeSubdomains();
	}

	size_t SPHSimulation::getParticleCount() const
	{
		return particles->getParticleCount();
	}

	void SPHSimulation::setGravity(float x, float y)
	{
		physics->setGravity(x, y);
	}

	void SPHSimulation::setViscosity(float v)
	{
		physics->setViscosity(v);
	}

	void SPHSimulation::setGasConstant(float k)
	{
		physics->setGasConstant(k);
	}

	void SPHSimulation::setRestDensity(float d)
	{
		physics->setRestDensity(d);
	}

	void SPHSimulation::applyMouseForce(const sf::Vector2f &mousePos, float strength)
	{
		// Radius of influence for mouse interaction
		const float radiusOfInfluence = 50.0f;
		const float radiusSqr = radiusOfInfluence * radiusOfInfluence;

		for (auto *particle : particles->getParticles())
		{
			sf::Vector2f particlePos = particle->getPosition();
			sf::Vector2f direction = particlePos - mousePos;
			float distanceSqr = direction.x * direction.x + direction.y * direction.y;

			// Only apply force if particle is within radius of influence
			if (distanceSqr < radiusSqr)
			{
				// Calculate force based on distance (stronger when closer)
				float distance = std::sqrt(distanceSqr);

				// Normalize direction and scale by strength and distance factor
				if (distance > 0.1f)
				{ // Prevent division by zero or very small values
					direction /= distance;
					float forceFactor = strength * (1.0f - distance / radiusOfInfluence);

					// Apply force directly to velocity for immediate response
					sf::Vector2f currentVel = particle->getVelocity();
					particle->setVelocity(currentVel + direction * forceFactor);
				}
			}
		}
	}

	void SPHSimulation::setLoadBalancingEnabled(bool enabled)
	{
		parallelExecutor->setLoadBalancingEnabled(enabled);
		std::cout << "Load balancing " << (enabled ? "enabled" : "disabled") << std::endl;
	}

	bool SPHSimulation::isLoadBalancingEnabled() const
	{
		return parallelExecutor->isLoadBalancingEnabled();
	}

	void SPHSimulation::setLoadBalanceThreshold(float threshold)
	{
		parallelExecutor->setLoadBalanceThreshold(threshold);
	}

	void SPHSimulation::setLoadBalanceInterval(int frames)
	{
		parallelExecutor->setLoadBalanceInterval(frames);
	}

	void SPHSimulation::setVisualizeLoadBalance(bool enabled)
	{
		renderer->setVisualizeLoadBalance(enabled);
		std::cout << "Load balance visualization " << (enabled ? "enabled" : "disabled") << std::endl;
	}

	bool SPHSimulation::isVisualizeLoadBalance() const
	{
		return renderer->isVisualizeLoadBalance();
	}

} // namespace sph