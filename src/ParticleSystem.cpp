#include "ParticleSystem.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace sph
{

	ParticleSystem::ParticleSystem(float width, float height, float smoothingRadius)
		: width(width),
		  height(height),
		  smoothingRadius(smoothingRadius)
	{
		// Create Grid for spatial partitioning
		grid = new Grid(width, height, smoothingRadius);
	}

	ParticleSystem::~ParticleSystem()
	{
		reset();
		delete grid;
	}

	Particle *ParticleSystem::addParticle(float x, float y)
	{
		Particle *particle = new Particle(x, y);
		particles.push_back(particle);
		return particle;
	}

	void ParticleSystem::addParticles(int count)
	{
		if (count <= 0)
			return;

		// Add particles in a small area near the center of the screen
		const float centerX = width / 2.0f;
		const float centerY = height / 3.0f; // Upper third of the screen
		const float radius = std::min(width, height) / 8.0f;

		for (int i = 0; i < count; ++i)
		{
			// Generate random position within a circle
			float angle = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 2.0f * M_PI;
			float distance = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * radius;

			float x = centerX + cos(angle) * distance;
			float y = centerY + sin(angle) * distance;

			// Ensure within bounds
			x = std::max(5.0f, std::min(width - 5.0f, x));
			y = std::max(5.0f, std::min(height - 5.0f, y));

			addParticle(x, y);
		}

		std::cout << "Added " << count << " particles, total: " << particles.size() << std::endl;
	}

	void ParticleSystem::removeParticles(int count)
	{
		if (count <= 0)
			return;

		// Remove the last N particles (or all if count > particles.size())
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

	void ParticleSystem::initialize(int count)
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

	void ParticleSystem::initializeDamBreak(int count, float damWidth, float damHeight)
	{
		// Clear any existing particles first
		reset();

		// Clamp dam dimensions to reasonable values
		damWidth = 0.1f;
		damHeight = 1.0f;

		// Calculate dam dimensions in world units
		float damWidthUnits = width * damWidth;
		float damHeightUnits = height * damHeight;

		// Calculate number of particles per dimension to achieve desired count
		// while maintaining reasonable density
		float particleSpacing = 10.0f;
		int columns = static_cast<int>(damWidthUnits / particleSpacing);
		int rows = static_cast<int>(damHeightUnits / particleSpacing);

		// Adjust particle spacing to fit exact count if possible
		if (columns * rows > count)
		{
			// Reduce density if too many particles
			float ratio = std::sqrt(static_cast<float>(count) / (columns * rows));
			columns = static_cast<int>(columns * ratio);
			rows = static_cast<int>(rows * ratio);
			std::cout << "Adjusted particle grid to " << columns << "x" << rows << " = "
					  << (columns * rows) << " particles (requested " << count << ")" << std::endl;
		}

		float spacingX = damWidthUnits / columns;
		float spacingY = damHeightUnits / rows;

		// Small jitter for more natural appearance
		float jitter = particleSpacing * 0.1f;

		std::cout << "Creating dam break scenario with " << columns << "x" << rows << " particles" << std::endl;

		// Create the dam block of particles (positioned in the left side of the domain)
		for (int y = 0; y < rows; ++y)
		{
			for (int x = 0; x < columns; ++x)
			{
				float posX = x * spacingX;
				// Invert y-position to start from bottom
				float posY = height - (y + 1) * spacingY;

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

		std::cout << "Dam break scenario created with " << particles.size() << " particles" << std::endl;
	}

	void ParticleSystem::reset()
	{
		for (auto *particle : particles)
		{
			delete particle;
		}
		particles.clear();
	}

	void ParticleSystem::updateGrid()
	{
		grid->clear();
		for (auto *particle : particles)
		{
			grid->insertParticle(particle);
		}
	}

} // namespace sph