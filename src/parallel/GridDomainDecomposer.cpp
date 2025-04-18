#include "parallel/GridDomainDecomposer.hpp"
#include <cmath>
#include <iostream>

namespace sph
{
	namespace parallel
	{
		std::vector<std::unique_ptr<Subdomain>> GridDomainDecomposer::createDecomposition(
			float width, float height, int numSubdomains)
		{
			std::vector<std::unique_ptr<Subdomain>> subdomains;

			// Calculate optimal grid dimensions
			int gridCols, gridRows;
			calculateGridDimensions(numSubdomains, width, height, gridCols, gridRows);

			// Calculate cell dimensions
			float cellWidth = width / gridCols;
			float cellHeight = height / gridRows;

			// Create subdomains in grid layout
			int id = 0;
			for (int row = 0; row < gridRows; row++)
			{
				for (int col = 0; col < gridCols; col++)
				{
					float subdomainX = col * cellWidth;
					float subdomainY = row * cellHeight;

					subdomains.push_back(
						std::make_unique<Subdomain>(
							subdomainX, subdomainY, cellWidth, cellHeight, id++));
				}
			}

			std::cout << "Created " << subdomains.size() << " subdomains in "
					  << gridCols << "x" << gridRows << " grid" << std::endl;

			return subdomains;
		}

		void GridDomainDecomposer::assignParticlesToSubdomains(
			const std::vector<Particle *> &particles,
			std::vector<std::unique_ptr<Subdomain>> &subdomains)
		{
			// Clear all subdomains first
			for (auto &subdomain : subdomains)
			{
				subdomain->clearParticles();
			}

			// Assign each particle to its corresponding subdomain
			for (auto *particle : particles)
			{
				sf::Vector2f pos = particle->getPosition();

				// Find subdomain that contains this particle
				for (auto &subdomain : subdomains)
				{
					if (subdomain->containsPoint(pos.x, pos.y))
					{
						subdomain->addParticle(particle);
						break;
					}
				}
			}
		}

		bool GridDomainDecomposer::updateDecomposition(
			std::vector<std::unique_ptr<Subdomain>> &subdomains,
			const std::vector<Particle *> &particles)
		{
			// Grid decomposition is static, so no updates needed
			// Just return false to indicate no changes were made
			return false;
		}

		void GridDomainDecomposer::calculateGridDimensions(
			int numSubdomains, float width, float height,
			int &gridCols, int &gridRows) const
		{
			// Calculate grid dimensions that maintain aspect ratio
			float aspectRatio = width / height;

			// Try to make cells as square as possible by considering aspect ratio
			gridCols = std::round(std::sqrt(numSubdomains * aspectRatio));
			gridRows = std::round(std::sqrt(numSubdomains / aspectRatio));

			// Ensure we get the requested number of subdomains
			while (gridCols * gridRows < numSubdomains)
			{
				if (gridCols / float(gridRows) < aspectRatio)
				{
					gridCols++;
				}
				else
				{
					gridRows++;
				}
			}

			// If we ended up with too many, reduce as needed
			while (gridCols * gridRows > numSubdomains && gridCols > 1 && gridRows > 1)
			{
				if (gridCols / float(gridRows) > aspectRatio)
				{
					gridCols--;
				}
				else
				{
					gridRows--;
				}
			}
		}
	}
}