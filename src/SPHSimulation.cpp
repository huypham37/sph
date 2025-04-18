#include "SPHSimulation.hpp"
#include <iostream>
#include <chrono>

namespace sph
{

	SPHSimulation::SPHSimulation(float width, float height)
		: width(width),
		  height(height),
		  smoothingRadius(16.0f)
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
		// Update grid for spatial partitioning
		particles->updateGrid();

		if (parallelExecutor->isParallelizationEnabled() && particles->getParticleCount() > 0)
		{
			// Update domain decomposition for parallel processing
			parallelExecutor->updateDecomposition(particles->getParticles());

			// Execute physics computations in parallel
			auto densityPressureTask = [this](const std::vector<Particle *> &localParticles,
											  const std::vector<Particle *> &ghostParticles)
			{
				// Process density and pressure using local and ghost particles
				for (auto *particle : localParticles)
				{
					float density = 0.0f;
					auto neighbors = particles->getGrid()->getNeighbors(particle, smoothingRadius);

					// Consider both regular neighbors and ghost particles
					auto processNeighbor = [&](Particle *neighbor)
					{
						sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
						float distSqr = r.x * r.x + r.y * r.y;

						if (distSqr < smoothingRadius * smoothingRadius)
						{
							density += neighbor->getMass() * physics->kernelPoly6(distSqr);
						}
					};

					// Process regular neighbors
					for (auto *neighbor : neighbors)
					{
						processNeighbor(neighbor);
					}

					// Process ghost neighbors
					for (auto *ghost : ghostParticles)
					{
						processNeighbor(ghost);
					}

					// Update particle density and pressure
					particle->setDensity(std::max(density, physics->getRestDensity()));
					float pressure = physics->getGasConstant() * (particle->getDensity() - physics->getRestDensity());
					particle->setPressure(std::max(0.0f, pressure));
				}
			};

			auto forcesTask = [this](const std::vector<Particle *> &localParticles,
									 const std::vector<Particle *> &ghostParticles)
			{
				// Process forces using local and ghost particles
				for (auto *particle : localParticles)
				{
					sf::Vector2f pressureForce = {0.0f, 0.0f};
					sf::Vector2f viscosityForce = {0.0f, 0.0f};

					auto neighbors = particles->getGrid()->getNeighbors(particle, smoothingRadius);

					auto processNeighbor = [&](Particle *neighbor)
					{
						// Skip self
						if (particle == neighbor)
							return;

						sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
						float dist = std::sqrt(r.x * r.x + r.y * r.y);

						if (dist > 0.0f && dist < smoothingRadius)
						{
							// Normalized direction
							sf::Vector2f dir = r / dist;

							// Pressure force
							float pressureTerm = particle->getPressure() / (particle->getDensity() * particle->getDensity()) +
												 neighbor->getPressure() / (neighbor->getDensity() * neighbor->getDensity());
							pressureForce -= neighbor->getMass() * pressureTerm * physics->kernelGradSpiky(dist, dir);

							// Viscosity force
							sf::Vector2f velocityDiff = neighbor->getVelocity() - particle->getVelocity();
							viscosityForce += physics->getViscosity() * neighbor->getMass() *
											  (velocityDiff / neighbor->getDensity()) *
											  physics->kernelViscosityLaplacian(dist);
						}
					};

					// Process regular neighbors
					for (auto *neighbor : neighbors)
					{
						processNeighbor(neighbor);
					}

					// Process ghost neighbors
					for (auto *ghost : ghostParticles)
					{
						processNeighbor(ghost);
					}

					// Gravity force
					sf::Vector2f gravityForce = physics->getGravity() * particle->getMass();

					// Total acceleration
					sf::Vector2f acceleration = (pressureForce + viscosityForce + gravityForce) / particle->getMass();
					particle->setAcceleration(acceleration);
				}
			};

			auto integrateTask = [dt](const std::vector<Particle *> &localParticles,
									  const std::vector<Particle *> &)
			{
				// Integrate particle positions and velocities
				for (auto *particle : localParticles)
				{
					// Semi-implicit Euler integration
					sf::Vector2f velocity = particle->getVelocity() + particle->getAcceleration() * dt;
					particle->setVelocity(velocity);
					particle->setPosition(particle->getPosition() + velocity * dt);
				}
			};

			// Execute tasks in parallel
			parallelExecutor->executeParallel(densityPressureTask, particles->getParticles());
			parallelExecutor->executeParallel(forcesTask, particles->getParticles());
			parallelExecutor->executeParallel(integrateTask, particles->getParticles());
		}
		else
		{
			// Sequential execution path
			physics->computeDensityPressure(particles->getParticles(), particles->getGrid());
			physics->computeForces(particles->getParticles(), particles->getGrid());
			physics->integrate(particles->getParticles(), dt);
		}

		// Collision resolution is always sequential
		physics->resolveCollisions(particles->getParticles(), particles->getGrid(), width, height);

		// Update particles for rendering
		for (auto *particle : particles->getParticles())
		{
			particle->update(dt);
		}
	}

	void SPHSimulation::draw(sf::RenderWindow &window)
	{
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
		particles->initialize(count);
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