#include "ParticleSystem.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include "SPHConfig.hpp"


namespace sph
{

	ParticleSystem::ParticleSystem(float width, float height, float smoothingRadius)
		: width(width),
		  height(height),
		  smoothingRadius(Config::SMOOTHING_RADIUS)
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



	void ParticleSystem::initializeDamBreak(int count)
	{
		// Clear any existing particles first
		reset();

		std::pair<int, int> ratio = Config::DAM_RATIO;
		float r = Config::PARTICLE_RADIUS;

		
		// Calculate dam dimensions in world units
		float damWidthUnits = ratio.second * 2 * r;
		float damHeightUnits = ratio.first * 2 * r;

		// Calculate number of particles per dimension to achieve desired count
		// while maintaining reasonable density
		float particleSpacing = Config::DP;
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