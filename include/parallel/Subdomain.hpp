#pragma once

#include <vector>
#include <memory>
#include <map>
#include "Particle.hpp"

namespace sph
{
	namespace parallel
	{
		/**
		 * @brief Enumeration of possible boundary regions in a 2D domain
		 */
		enum class BoundaryRegion
		{
			LEFT,
			RIGHT,
			TOP,
			BOTTOM,
			TOP_LEFT,
			TOP_RIGHT,
			BOTTOM_LEFT,
			BOTTOM_RIGHT
		};

		/**
		 * @brief Represents a subdomain in the simulation space
		 *
		 * A Subdomain is a rectangular section of the simulation domain that contains
		 * a subset of particles and is processed by a single thread/process.
		 */
		class Subdomain
		{
		public:
			Subdomain(float x, float y, float width, float height, int id);
			~Subdomain() = default;

			// Geometry access
			float getX() const { return x; }
			float getY() const { return y; }
			float getWidth() const { return width; }
			float getHeight() const { return height; }
			int getId() const { return id; }

			// Particle management
			void addParticle(Particle *particle);
			void removeParticle(Particle *particle);
			void clearParticles();
			const std::vector<Particle *> &getParticles() const { return particles; }

			// Ghost particles (from neighboring subdomains)
			void addGhostParticle(Particle *particle);
			void clearGhostParticles();
			const std::vector<Particle *> &getGhostParticles() const { return ghostParticles; }

			// Boundary region management
			void updateBoundaryRegions(float ghostRegionWidth);
			const std::vector<Particle *> &getParticlesInBoundaryRegion(BoundaryRegion region) const;
			std::vector<BoundaryRegion> getSharedBoundaryRegions(const Subdomain &other, float ghostRegionWidth) const;

			// Performance metrics
			float getLastComputationTime() const { return lastComputationTime; }
			void setLastComputationTime(float time) { lastComputationTime = time; }

			// Check if a point is within this subdomain
			bool containsPoint(float px, float py) const;

		private:
			float x;	  // Left coordinate of subdomain
			float y;	  // Top coordinate of subdomain
			float width;  // Width of subdomain
			float height; // Height of subdomain
			int id;		  // Unique identifier

			std::vector<Particle *> particles;									 // Particles owned by this subdomain
			std::vector<Particle *> ghostParticles;								 // Particles from neighboring subdomains
			std::map<BoundaryRegion, std::vector<Particle *>> boundaryParticles; // Particles near boundaries

			float lastComputationTime; // Time taken to process this subdomain (for load balancing)
		};

	} // namespace parallel
} // namespace sph