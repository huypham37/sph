#include "SPHSolver.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <omp.h> // Added OpenMP header

SPHSolver::SPHSolver(float width, float height)
	: width(width),
	  height(height),
	  gravity({0.0f, 9.8f}), // Using braced initialization for SFML 3.0
	  h(16.0f),
	  h2(h * h),
	  viscosityCoefficient(0.1f),
	  gasConstant(200.0f),
	  restDensity(1000.0f),
	  boundaryDamping(0.5f)
{
	grid = new Grid(width, height, h);
}

SPHSolver::~SPHSolver()
{
	reset();
	delete grid;
}

void SPHSolver::update(float dt)
{
	// Reset grid for spatial partitioning
	grid->clear();

	// Insert particles into grid
	for (auto *particle : particles)
	{
		grid->insertParticle(particle);
	}

	// Calculate density and pressure for each particle
	computeDensityPressure();

	// Calculate forces (pressure, viscosity, gravity)
	computeForces();

	// Integrate motion (update positions)
	integrate(dt);

	// Handle boundary collisions
	resolveCollisions();

	// Update particles for rendering
	for (auto *particle : particles)
	{
		particle->update(dt);
	}
}

void SPHSolver::draw(sf::RenderWindow &window)
{
	for (auto *particle : particles)
	{
		particle->draw(window);
	}
}

void SPHSolver::addParticle(float x, float y)
{
	particles.push_back(new Particle(x, y));
}

void SPHSolver::reset()
{
	for (auto *particle : particles)
	{
		delete particle;
	}
	particles.clear();
}

void SPHSolver::computeDensityPressure()
{
	for (auto *particle : particles)
	{
		// Get neighbors within smoothing radius h
		auto neighbors = grid->getNeighbors(particle, h);

		// Calculate density
		float density = 0.0f;
		for (auto *neighbor : neighbors)
		{
			sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
			float distSqr = r.x * r.x + r.y * r.y;

			if (distSqr < h2)
			{
				// Apply smoothing kernel (Poly6)
				density += neighbor->getMass() * kernelPoly6(distSqr);
			}
		}

		// Update particle density
		particle->setDensity(std::max(density, restDensity));

		// Calculate pressure using equation of state (ideal gas)
		float pressure = gasConstant * (particle->getDensity() - restDensity);
		particle->setPressure(std::max(0.0f, pressure));
	}
}

void SPHSolver::computeForces()
{
	for (auto *particle : particles)
	{
		sf::Vector2f pressureForce = {0.0f, 0.0f};	// Using braced initialization
		sf::Vector2f viscosityForce = {0.0f, 0.0f}; // Using braced initialization

		auto neighbors = grid->getNeighbors(particle, h);

		for (auto *neighbor : neighbors)
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

				// Pressure force: -∇p = -grad(p) = -mass * (p_i/rho_i^2 + p_j/rho_j^2) * grad_W
				float pressureTerm = particle->getPressure() / (particle->getDensity() * particle->getDensity()) +
									 neighbor->getPressure() / (neighbor->getDensity() * neighbor->getDensity());

				pressureForce -= neighbor->getMass() * pressureTerm * kernelGradSpiky(dist, dir);

				// Viscosity force: μ∇^2v = μ * mass * (v_j - v_i)/rho_j * laplacian_W
				sf::Vector2f velocityDiff = neighbor->getVelocity() - particle->getVelocity();
				viscosityForce += viscosityCoefficient * neighbor->getMass() *
								  (velocityDiff / neighbor->getDensity()) *
								  kernelViscosityLaplacian(dist);
			}
		}

		// Gravity force
		sf::Vector2f gravityForce = gravity * particle->getMass();

		// Calculate total acceleration: F/mass
		sf::Vector2f acceleration = (pressureForce + viscosityForce + gravityForce) / particle->getMass();

		// Update particle acceleration
		particle->setAcceleration(acceleration);
	}
}

void SPHSolver::integrate(float dt)
{
// Add debug output to verify OpenMP is working
#pragma omp parallel
	{
#pragma omp single
		{
			std::cout << "OpenMP is using " << omp_get_num_threads() << " threads" << std::endl;
		}
	}

#pragma omp parallel for
	for (size_t i = 0; i < particles.size(); i++)
	{
		Particle *particle = particles[i];
		// Semi-implicit Euler integration
		sf::Vector2f velocity = particle->getVelocity() + particle->getAcceleration() * dt;
		particle->setVelocity(velocity);
		particle->setPosition(particle->getPosition() + velocity * dt);
	}
}

void SPHSolver::resolveCollisions()
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
float SPHSolver::kernelPoly6(float distSqr)
{
	if (distSqr >= h2)
		return 0.0f;

	const float coeff = 315.0f / (64.0f * M_PI * std::pow(h, 9));
	float h2_r2 = h2 - distSqr;
	return coeff * h2_r2 * h2_r2 * h2_r2;
}

sf::Vector2f SPHSolver::kernelGradSpiky(float dist, const sf::Vector2f &dir)
{
	if (dist >= h || dist <= 0.0f)
		return sf::Vector2f{0.0f, 0.0f}; // Using braced initialization for SFML 3.0

	const float coeff = -45.0f / (M_PI * std::pow(h, 6));
	float h_r = h - dist;
	return coeff * h_r * h_r * dir;
}

float SPHSolver::kernelViscosityLaplacian(float dist)
{
	if (dist >= h)
		return 0.0f;

	const float coeff = 45.0f / (M_PI * std::pow(h, 6));
	return coeff * (h - dist);
}

void SPHSolver::initializeDefaultParticles(int count)
{
	// Clear any existing particles first
	reset();

	// Calculate grid parameters for even distribution
	int columns = static_cast<int>(sqrt(count * width / height));
	int rows = static_cast<int>(count / columns) + 1;

	float spacingX = width / (columns + 1);
	float spacingY = height / (rows + 1) * 0.5f; // Use only top half of the screen

	// Create a grid of particles
	for (int y = 0; y < rows; ++y)
	{
		for (int x = 0; x < columns; ++x)
		{
			// Add some randomness to positions for more natural look
			float offsetX = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 5.0f - 2.5f;
			float offsetY = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 5.0f - 2.5f;

			float posX = (x + 1) * spacingX + offsetX;
			float posY = (y + 1) * spacingY + offsetY;

			// Ensure we don't exceed the desired particle count
			if (particles.size() < static_cast<size_t>(count))
			{
				addParticle(posX, posY);
			}
			else
			{
				return;
			}
		}
	}
}

void SPHSolver::applyMouseForce(const sf::Vector2f &mousePos, float strength)
{
	// Radius of influence for mouse interaction
	const float radiusOfInfluence = 50.0f;
	const float radiusSqr = radiusOfInfluence * radiusOfInfluence;

	for (auto *particle : particles)
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