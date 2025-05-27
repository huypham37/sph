#include "SPHSimulation.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <cmath> // Add cmath header for std::sqrt
#ifdef _OPENMP
#include <omp.h> // Include OpenMP header for thread functions
#endif

namespace sph
{

	SPHSimulation::SPHSimulation(float width, float height)
		: width(width),
		  height(height),
		  smoothingRadius(Config::SMOOTHING_RADIUS),
		  framesSinceLastMetricUpdate(0),
		  particlesPerSecond(0.0),
		  timeStepsPerSecond(0.0),
		  simToPhysicalRatio(0.0),
		  avgTimePerStep(0.0),
		  totalSimulationTime(0.0),
		  totalPhysicalTime(0.0)
	{
		// Create component instances
		particles = std::make_unique<ParticleSystem>(width, height, smoothingRadius);
		physics = std::make_unique<SPHPhysics>();
		physics->setSmoothingRadius(smoothingRadius);
#ifndef HEADLESS_MODE
		renderer = std::make_unique<Renderer>();
#endif

		lastMetricUpdateTime = std::chrono::high_resolution_clock::now();

		std::cout << "SPH Simulation initialized with dimensions " << width << "x" << height << std::endl;
	}

	SPHSimulation::~SPHSimulation()
	{
		// Unique pointers will clean up automatically
	}

	void SPHSimulation::update(float dt)
	{
		// Start timing this step
		auto stepStartTime = std::chrono::high_resolution_clock::now();

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
			physics->resolveCollisionsForParticle(particle, width, height);
		}

		// Calculate metrics
		auto stepEndTime = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> stepDuration = stepEndTime - stepStartTime;
		double stepTime = stepDuration.count();

		// Update metrics
		framesSinceLastMetricUpdate++;
		totalSimulationTime += stepTime;
		totalPhysicalTime += dt;
		simToPhysicalRatio = totalSimulationTime / totalPhysicalTime;

		stepTimes.push_back(stepTime);
		if (stepTimes.size() > 100)
			stepTimes.erase(stepTimes.begin()); // Keep last 100 measurements

		// Calculate average time per step
		avgTimePerStep = 0;
		for (double time : stepTimes)
			avgTimePerStep += time;
		avgTimePerStep /= stepTimes.size();

		// Update metrics every second
		auto now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = now - lastMetricUpdateTime;
		if (elapsed.count() >= 1.0)
		{ // Update metrics every second
			particlesPerSecond = particles->getParticleCount() * framesSinceLastMetricUpdate / elapsed.count();
			timeStepsPerSecond = framesSinceLastMetricUpdate / elapsed.count();

			framesSinceLastMetricUpdate = 0;
			lastMetricUpdateTime = now;
		}
	}

#ifndef HEADLESS_MODE
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
#endif

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

	void SPHSimulation::printPerformanceMetrics() const
	{
		std::cout << "===== Performance Metrics =====" << std::endl;
		std::cout << "Particles: " << particles->getParticleCount() << std::endl;
#ifdef _OPENMP
		std::cout << "OpenMP threads: " << omp_get_max_threads() << std::endl;
#else
		std::cout << "OpenMP: Not enabled" << std::endl;
#endif
		std::cout << "Particles/second: " << std::fixed << std::setprecision(1) << particlesPerSecond << std::endl;
		std::cout << "Timesteps/second: " << std::fixed << std::setprecision(2) << timeStepsPerSecond << std::endl;
		std::cout << "Timesteps/minute: " << std::fixed << std::setprecision(1) << timeStepsPerSecond * 60 << std::endl;
		std::cout << "Timesteps/hour: " << std::fixed << std::setprecision(1) << timeStepsPerSecond * 3600 << std::endl;
		std::cout << "Simulation:Physical time ratio: " << std::fixed << std::setprecision(3) << simToPhysicalRatio << "x" << std::endl;
		std::cout << "Avg time per step: " << std::fixed << std::setprecision(5) << avgTimePerStep * 1000 << " ms" << std::endl;
		std::cout << "=============================" << std::endl;
	}

} // namespace sph