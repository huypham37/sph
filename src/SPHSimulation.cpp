#include "SPHSimulation.hpp"
#include <iostream>
#include <chrono>

namespace sph
{

	SPHSimulation::SPHSimulation(float width, float height)
		: width(width),
		  height(height),
		  smoothingRadius(Config::SMOOTHING_RADIUS)
	{
		// Create component instances
		particles = std::make_unique<ParticleSystem>(width, height, smoothingRadius);
		physics = std::make_unique<SPHPhysics>();
		physics->setSmoothingRadius(smoothingRadius);
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

		// Update grid (needs to be sequential due to data structure updates)
		particles->updateGrid();

		// Substep loop for physics
		for (int step = 0; step < sub_steps; ++step)
		{
// Single OpenMP parallel region for the entire simulation step
#pragma omp parallel
			{
// Cache neighbors in parallel
#pragma omp for schedule(dynamic, 64)
				for (size_t i = 0; i < particles->getParticles().size(); ++i)
				{
					auto *particle = particles->getParticles()[i];
					particle->cachedNeighbors = particles->getGrid()->getNeighbors(particle, smoothingRadius * 2);
				}

// Density and pressure calculation phase
#pragma omp barrier
#pragma omp for schedule(dynamic, 64)
				for (size_t i = 0; i < particles->getParticles().size(); ++i)
				{
					auto *particle = particles->getParticles()[i];
					physics->computeDensityPressureForParticle(particle);
				}

// Force calculation phase
#pragma omp barrier
#pragma omp for schedule(dynamic, 100)
				for (size_t i = 0; i < particles->getParticles().size(); ++i)
				{
					auto *particle = particles->getParticles()[i];
					physics->computeForcesForParticle(particle);
				}

// Boundary forces
#pragma omp barrier
#pragma omp for schedule(dynamic, 100)
				for (size_t i = 0; i < particles->getParticles().size(); ++i)
				{
					auto *particle = particles->getParticles()[i];
					physics->computeBoundaryForcesForParticle(particle, width, height);
				}

// Integration
#pragma omp barrier
#pragma omp for schedule(dynamic, 100)
				for (size_t i = 0; i < particles->getParticles().size(); ++i)
				{
					auto *particle = particles->getParticles()[i];
					physics->integrateParticle(particle, sub_dt);
					particle->update(sub_dt);
				}
			} // End of parallel region
		}

// Collision resolution (can be parallel but separate from main computation)
#pragma omp parallel for schedule(dynamic, 100)
		for (size_t i = 0; i < particles->getParticles().size(); ++i)
		{
			auto *particle = particles->getParticles()[i];
			physics->resolveCollisionsForParticle(particle,width, height);
		}
	}

	void SPHSimulation::draw(sf::RenderWindow &window)
	{
		// Only update visuals for particles that are visible
		for (auto *particle : particles->getParticles())
		{
			particle->updateVisuals();
		}

		// Draw particles
		renderer->drawParticles(particles->getParticles(), window);
	}

	void SPHSimulation::reset()
	{
		particles->reset();
	}

	void SPHSimulation::initializeDefaultParticles(int count)
	{
		particles->initializeDamBreak(count);
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

} // namespace sph