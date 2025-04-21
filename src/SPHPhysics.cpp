#include "SPHPhysics.hpp"
#include "Grid.hpp"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <unordered_set>

namespace sph
{

	SPHPhysics::SPHPhysics()
		: h(0.2f),
		  h2(h * h),
		  viscosityCoefficient(0.1f),
		  gasConstant(200.0f),
		  restDensity(1000.0f),
		  boundaryDamping(0.5f)
	{
	}

	void SPHPhysics::setSmoothingRadius(float smoothingRadius)
	{
		h = smoothingRadius;
		h2 = h * h; // Pre-compute h^2 for optimization
	}

	void SPHPhysics::computeDensityPressure(const std::vector<Particle *> &particles, Grid *grid)
	{
		for (auto *particle : particles)
		{
			float density = 0.0f;

			// Use cached neighbors
			for (auto *neighbor : particle->cachedNeighbors)
			{
				sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
				float distSqr = r.x * r.x + r.y * r.y;

				if (distSqr < h * h)
				{
					density += neighbor->getMass() * kernelPoly6(distSqr);
				}
			}

			// Update particle density and pressure
			particle->setDensity(std::max(density, restDensity));
			float pressure = gasConstant * (particle->getDensity() - restDensity);
			particle->setPressure(std::max(0.0f, pressure));
		}
	}

	void SPHPhysics::computeForces(const std::vector<Particle *> &particles, Grid *grid)
	{
		for (auto *particle : particles)
		{
			sf::Vector2f pressureForce = {0.0f, 0.0f};
			sf::Vector2f viscosityForce = {0.0f, 0.0f};

			// Use cached neighbors
			for (auto *neighbor : particle->cachedNeighbors)
			{
				// Skip self
				if (particle == neighbor)
					continue;

				sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
				float dist = std::sqrt(r.x * r.x + r.y * r.y);

				if (dist > 0.0f && dist < h)
				{
					// Normalized direction
					sf::Vector2f dir = r / dist;

					// Pressure force
					float pressureTerm = particle->getPressure() / (particle->getDensity() * particle->getDensity()) +
										 neighbor->getPressure() / (neighbor->getDensity() * neighbor->getDensity());
					pressureForce -= neighbor->getMass() * pressureTerm * kernelGradSpiky(dist, dir);

					// Viscosity force
					sf::Vector2f velocityDiff = neighbor->getVelocity() - particle->getVelocity();
					viscosityForce += viscosityCoefficient * neighbor->getMass() *
									  (velocityDiff / neighbor->getDensity()) *
									  kernelViscosityLaplacian(dist);
				}
			}

			// Gravity force
			sf::Vector2f gravityForce = gravity * particle->getMass();

			// Total acceleration
			sf::Vector2f acceleration = (pressureForce + viscosityForce + gravityForce) / particle->getMass();
			particle->setAcceleration(acceleration);
		}
	}

	void SPHPhysics::integrate(const std::vector<Particle *> &particles, float dt)
	{
		for (auto *particle : particles)
		{
			// Semi-implicit Euler integration
			sf::Vector2f velocity = particle->getVelocity() + particle->getAcceleration() * dt;
			particle->setVelocity(velocity);
			particle->setPosition(particle->getPosition() + velocity * dt);
		}
	}

	void SPHPhysics::resolveCollisions(const std::vector<Particle *> &particles, Grid *grid, float width, float height)
	{
		constexpr float PARTICLE_RADIUS = 4.0f;	  // Should match the radius defined in Particle class
		constexpr float COLLISION_DAMPING = 0.5f; // Damping factor for particle collisions

		// First handle boundary collisions
		for (auto *particle : particles)
		{
			sf::Vector2f pos = particle->getPosition();
			sf::Vector2f vel = particle->getVelocity();

			// Boundary checks with damping
			// Left boundary
			if (pos.x < PARTICLE_RADIUS)
			{
				pos.x = PARTICLE_RADIUS;
				vel.x = -vel.x * boundaryDamping;
			}
			// Right boundary
			if (pos.x > width - PARTICLE_RADIUS)
			{
				pos.x = width - PARTICLE_RADIUS;
				vel.x = -vel.x * boundaryDamping;
			}
			// Top boundary
			if (pos.y < PARTICLE_RADIUS)
			{
				pos.y = PARTICLE_RADIUS;
				vel.y = -vel.y * boundaryDamping;
			}
			// Bottom boundary
			if (pos.y > height - PARTICLE_RADIUS)
			{
				pos.y = height - PARTICLE_RADIUS;
				vel.y = -vel.y * boundaryDamping;
			}

			particle->setPosition(pos);
			particle->setVelocity(vel);
		}

		// Then handle particle-particle collisions
		const float minDist = PARTICLE_RADIUS * 2.0f; // Minimum distance between particles
		const float minDistSq = minDist * minDist;

		// We'll use our grid to find potential collision pairs efficiently
		grid->clear();
		for (auto *particle : particles)
		{
			grid->insertParticle(particle);
		}

		// Check each particle against potential neighbors
		std::vector<std::tuple<Particle *, sf::Vector2f, sf::Vector2f>> updates;

		for (size_t i = 0; i < particles.size(); ++i)
		{
			Particle *p1 = particles[i];
			sf::Vector2f pos1 = p1->getPosition();

			// Get potential collision candidates from grid
			std::vector<Particle *> neighbors = grid->getNeighbors(p1, minDist * 1.5f); // slightly larger radius for safety

			for (auto *p2 : neighbors)
			{
				// Skip self-collision
				if (p1 == p2)
					continue;

				sf::Vector2f pos2 = p2->getPosition();
				sf::Vector2f delta = pos1 - pos2;
				float distSq = delta.x * delta.x + delta.y * delta.y;

				// If particles are overlapping
				if (distSq < minDistSq && distSq > 0.0001f)
				{ // Avoid division by zero
					float dist = std::sqrt(distSq);
					float penetration = minDist - dist;

					// Normalized collision normal
					sf::Vector2f normal = delta / dist;

					// Resolve positions - push both particles apart
					sf::Vector2f correction = normal * (penetration * 0.5f);

					// Velocity response - elastic collision
					sf::Vector2f v1 = p1->getVelocity();
					sf::Vector2f v2 = p2->getVelocity();

					// Relative velocity
					sf::Vector2f relVel = v1 - v2;
					float velAlongNormal = relVel.x * normal.x + relVel.y * normal.y;

					// Only separate if particles are moving toward each other
					if (velAlongNormal < 0)
					{
						// Apply impulse
						sf::Vector2f impulse = normal * (-velAlongNormal * COLLISION_DAMPING);

						updates.push_back(std::make_tuple(p1, pos1 + correction, v1 + impulse));
						updates.push_back(std::make_tuple(p2, pos2 - correction, v2 - impulse));
					}
					else
					{
						// Just correct positions if not colliding
						updates.push_back(std::make_tuple(p1, pos1 + correction, v1));
						updates.push_back(std::make_tuple(p2, pos2 - correction, v2));
					}
				}
			}
		}

		// Apply all position and velocity updates
		for (const auto &update : updates)
		{
			Particle *p = std::get<0>(update);
			p->setPosition(std::get<1>(update));
			p->setVelocity(std::get<2>(update));
		}
	}

	// SPH Kernel functions
	float SPHPhysics::kernelPoly6(float distSquared)
	{
		if (distSquared >= h2)
			return 0.0f;

		const float coeff = 315.0f / (64.0f * M_PI * std::pow(h, 9));
		float h2_r2 = h2 - distSquared;
		return coeff * h2_r2 * h2_r2 * h2_r2;
	}

	sf::Vector2f SPHPhysics::kernelGradSpiky(float dist, const sf::Vector2f &dir)
	{
		if (dist >= h || dist <= 0.0f)
			return {0.0f, 0.0f};

		const float coeff = -45.0f / (M_PI * std::pow(h, 6));
		float h_r = h - dist;
		return coeff * h_r * h_r * dir;
	}

	float SPHPhysics::kernelViscosityLaplacian(float dist)
	{
		if (dist >= h)
			return 0.0f;

		const float coeff = 45.0f / (M_PI * std::pow(h, 6));
		return coeff * (h - dist);
	}

} // namespace sph