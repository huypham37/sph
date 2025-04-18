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

		bool Subdomain::containsPoint(float px, float py) const
		{
			return (px >= x && px < x + width && py >= y && py < y + height);
		}
	}
}