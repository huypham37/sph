#include "SPHSolver.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

SPHSolver::SPHSolver(float width, float height)
	: width(width),
	  height(height),
	  physicsEngine(), // Initialize physics engine
	  fontLoaded(false)
{
	// Set physics parameters
	float initialSmoothingRadius = 16.0f;
	physicsEngine.setSmoothingRadius(initialSmoothingRadius);
	physicsEngine.setViscosity(0.1f);
	physicsEngine.setGasConstant(200.0f);
	physicsEngine.setRestDensity(1000.0f);
	physicsEngine.setBoundaryDamping(0.5f);
	physicsEngine.setGravity(0.0f, 9.8f); // Fixed to use correct signature

	// Initialize grid for spatial partitioning
	grid = new Grid(width, height, initialSmoothingRadius);

	// Try to load a system font for visualization if needed
	fontLoaded = font.openFromFile("/System/Library/Fonts/Helvetica.ttc");
	if (!fontLoaded)
	{
		std::cout << "Warning: Failed to load font for visualization" << std::endl;
	}
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

	// Sequential SPH computation using physicsEngine
	physicsEngine.computeDensityPressure(particles, grid);
	physicsEngine.computeForces(particles, grid);
	physicsEngine.integrate(particles, dt);
	physicsEngine.resolveCollisions(particles, grid, width, height);

	// Update particles for rendering (visual state update)
	for (auto *particle : particles)
	{
		particle->update(dt);
	}
}

void SPHSolver::draw(sf::RenderWindow &window)
{
	// Draw all particles
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

void SPHSolver::setSmoothingRadius(float radius)
{
	physicsEngine.setSmoothingRadius(radius);

	// Since Grid doesn't have updateCellSize method, we need to recreate the grid
	delete grid;
	grid = new Grid(width, height, radius);

	// Re-insert particles into the new grid
	for (auto *particle : particles)
	{
		grid->insertParticle(particle);
	}
}

void SPHSolver::initializeDefaultParticles(int count)
{
	reset();

	int columns = static_cast<int>(sqrt(count * width / height));
	int rows = static_cast<int>(count / columns) + 1;

	float spacingX = width / (columns + 1);
	float spacingY = height / (rows + 1) * 0.5f;

	for (int y = 0; y < rows; ++y)
	{
		for (int x = 0; x < columns; ++x)
		{
			float offsetX = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 5.0f - 2.5f;
			float offsetY = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 5.0f - 2.5f;

			float posX = (x + 1) * spacingX + offsetX;
			float posY = (y + 1) * spacingY + offsetY;

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
	const float radiusOfInfluence = 50.0f;
	const float radiusSqr = radiusOfInfluence * radiusOfInfluence;

	for (auto *particle : particles)
	{
		sf::Vector2f particlePos = particle->getPosition();
		sf::Vector2f direction = particlePos - mousePos;
		float distanceSqr = direction.x * direction.x + direction.y * direction.y;

		if (distanceSqr < radiusSqr)
		{
			float distance = std::sqrt(distanceSqr);

			if (distance > 0.1f)
			{
				direction /= distance;
				float forceFactor = strength * (1.0f - distance / radiusOfInfluence);

				sf::Vector2f currentVel = particle->getVelocity();
				particle->setVelocity(currentVel + direction * forceFactor);
			}
		}
	}
}

void SPHSolver::addParticles(int count)
{
	if (count <= 0)
		return;

	const float centerX = width / 2.0f;
	const float centerY = height / 3.0f;
	const float radius = std::min(width, height) / 8.0f;

	for (int i = 0; i < count; ++i)
	{
		float angle = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 2.0f * M_PI;
		float distance = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * radius;

		float x = centerX + cos(angle) * distance;
		float y = centerY + sin(angle) * distance;

		x = std::max(5.0f, std::min(width - 5.0f, x));
		y = std::max(5.0f, std::min(height - 5.0f, y));

		addParticle(x, y);
	}

	std::cout << "Added " << count << " particles, total: " << particles.size() << std::endl;
}

void SPHSolver::removeParticles(int count)
{
	if (count <= 0)
		return;

	count = std::min(count, static_cast<int>(particles.size()));

	for (int i = 0; i < count; ++i)
	{
		if (!particles.empty())
		{
			delete particles.back();
			particles.pop_back();
		}
	}

	std::cout << "Removed " << count << " particles, total: " << particles.size() << std::endl;
}