#include "parallel/Subdomain.hpp"
#include <algorithm>

namespace sph
{
	namespace parallel
	{
		Subdomain::Subdomain(float x, float y, float width, float height, int id)
			: x(x), y(y), width(width), height(height), id(id), lastComputationTime(0.0f)
		{
		}

		void Subdomain::addParticle(Particle *particle)
		{
			// Check if particle is already in the subdomain
			auto it = std::find(particles.begin(), particles.end(), particle);
			if (it == particles.end())
			{
				particles.push_back(particle);
			}
		}

		void Subdomain::removeParticle(Particle *particle)
		{
			auto it = std::find(particles.begin(), particles.end(), particle);
			if (it != particles.end())
			{
				particles.erase(it);
			}
		}

		void Subdomain::clearParticles()
		{
			particles.clear();
		}

		void Subdomain::addGhostParticle(Particle *particle)
		{
			// Check if particle is already in ghost particles
			auto it = std::find(ghostParticles.begin(), ghostParticles.end(), particle);
			if (it == ghostParticles.end())
			{
				ghostParticles.push_back(particle);
			}
		}

		void Subdomain::clearGhostParticles()
		{
			ghostParticles.clear();
		}

		void Subdomain::updateBoundaryRegions(float ghostRegionWidth)
		{
			// Clear previous boundary particles
			boundaryParticles.clear();

			// Calculate boundary region extents
			float leftBound = x;
			float rightBound = x + width;
			float topBound = y;
			float bottomBound = y + height;

			// Initialize boundary region containers
			for (int i = 0; i < 8; i++)
			{
				boundaryParticles[static_cast<BoundaryRegion>(i)] = std::vector<Particle *>();
			}

			// Categorize each particle into appropriate boundary regions
			for (auto *particle : particles)
			{
				sf::Vector2f pos = particle->getPosition();

				// Check if particle is in each boundary region
				bool inLeftRegion = (pos.x >= leftBound && pos.x < leftBound + ghostRegionWidth);
				bool inRightRegion = (pos.x >= rightBound - ghostRegionWidth && pos.x < rightBound);
				bool inTopRegion = (pos.y >= topBound && pos.y < topBound + ghostRegionWidth);
				bool inBottomRegion = (pos.y >= bottomBound - ghostRegionWidth && pos.y < bottomBound);

				// Add to appropriate regions
				if (inLeftRegion)
				{
					boundaryParticles[BoundaryRegion::LEFT].push_back(particle);

					if (inTopRegion)
					{
						boundaryParticles[BoundaryRegion::TOP_LEFT].push_back(particle);
					}
					if (inBottomRegion)
					{
						boundaryParticles[BoundaryRegion::BOTTOM_LEFT].push_back(particle);
					}
				}

				if (inRightRegion)
				{
					boundaryParticles[BoundaryRegion::RIGHT].push_back(particle);

					if (inTopRegion)
					{
						boundaryParticles[BoundaryRegion::TOP_RIGHT].push_back(particle);
					}
					if (inBottomRegion)
					{
						boundaryParticles[BoundaryRegion::BOTTOM_RIGHT].push_back(particle);
					}
				}

				if (inTopRegion && !inLeftRegion && !inRightRegion)
				{
					boundaryParticles[BoundaryRegion::TOP].push_back(particle);
				}

				if (inBottomRegion && !inLeftRegion && !inRightRegion)
				{
					boundaryParticles[BoundaryRegion::BOTTOM].push_back(particle);
				}
			}
		}

		const std::vector<Particle *> &Subdomain::getParticlesInBoundaryRegion(BoundaryRegion region) const
		{
			auto it = boundaryParticles.find(region);
			if (it != boundaryParticles.end())
			{
				return it->second;
			}

			static std::vector<Particle *> emptyVector;
			return emptyVector;
		}

		std::vector<BoundaryRegion> Subdomain::getSharedBoundaryRegions(const Subdomain &other, float ghostRegionWidth) const
		{
			std::vector<BoundaryRegion> sharedRegions;

			// Calculate extended boundaries
			float thisLeft = x - ghostRegionWidth;
			float thisRight = x + width + ghostRegionWidth;
			float thisTop = y - ghostRegionWidth;
			float thisBottom = y + height + ghostRegionWidth;

			float otherLeft = other.getX();
			float otherRight = other.getX() + other.getWidth();
			float otherTop = other.getY();
			float otherBottom = other.getY() + other.getHeight();

			// Check for overlaps
			bool otherOnLeft = (otherRight > thisLeft && otherLeft < x);
			bool otherOnRight = (otherLeft < thisRight && otherRight > x + width);
			bool otherOnTop = (otherBottom > thisTop && otherTop < y);
			bool otherOnBottom = (otherTop < thisBottom && otherBottom > y + height);

			// Add appropriate regions
			if (otherOnLeft)
			{
				sharedRegions.push_back(BoundaryRegion::LEFT);
				if (otherOnTop)
					sharedRegions.push_back(BoundaryRegion::TOP_LEFT);
				if (otherOnBottom)
					sharedRegions.push_back(BoundaryRegion::BOTTOM_LEFT);
			}

			if (otherOnRight)
			{
				sharedRegions.push_back(BoundaryRegion::RIGHT);
				if (otherOnTop)
					sharedRegions.push_back(BoundaryRegion::TOP_RIGHT);
				if (otherOnBottom)
					sharedRegions.push_back(BoundaryRegion::BOTTOM_RIGHT);
			}

			if (otherOnTop && !otherOnLeft && !otherOnRight)
			{
				sharedRegions.push_back(BoundaryRegion::TOP);
			}

			if (otherOnBottom && !otherOnLeft && !otherOnRight)
			{
				sharedRegions.push_back(BoundaryRegion::BOTTOM);
			}

			return sharedRegions;
		}

		bool Subdomain::containsPoint(float px, float py) const
		{
			return (px >= x && px < x + width && py >= y && py < y + height);
		}
	}
}