#include "SPHPhysics.hpp"
#include "Grid.hpp"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <unordered_set>
#include <omp.h>
#include <iostream>

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
        h2 = h * h;
    }

    void SPHPhysics::computeDensityPressure(const std::vector<Particle *> &particles, Grid *grid)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < particles.size(); ++i)
        {
            float density = 0.0f;
            auto *particle = particles[i];

            for (auto *neighbor : particle->cachedNeighbors)
            {
                sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
                float distSqr = r.x * r.x + r.y * r.y;

                if (distSqr < h * h)
                {
                    density += neighbor->getMass() * kernelPoly6(distSqr);
                }
            }

            particle->setDensity(std::max(density, restDensity));
            float pressure = gasConstant * (particle->getDensity() - restDensity);
            particle->setPressure(std::max(0.0f, pressure));
        }
    }

    void SPHPhysics::computeForces(const std::vector<Particle *> &particles, Grid *grid)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < particles.size(); ++i)
        {
            auto *particle = particles[i];
            sf::Vector2f pressureForce = {0.0f, 0.0f};
            sf::Vector2f viscosityForce = {0.0f, 0.0f};

            for (auto *neighbor : particle->cachedNeighbors)
            {
                if (particle == neighbor)
                    continue;

                sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
                float dist = std::sqrt(r.x * r.x + r.y * r.y);

                if (dist > 0.0f && dist < h)
                {
                    sf::Vector2f dir = r / dist;

                    float pressureTerm = particle->getPressure() / (particle->getDensity() * particle->getDensity()) +
                                         neighbor->getPressure() / (neighbor->getDensity() * neighbor->getDensity());
                    pressureForce -= neighbor->getMass() * pressureTerm * kernelGradSpiky(dist, dir);

                    sf::Vector2f velocityDiff = neighbor->getVelocity() - particle->getVelocity();
                    viscosityForce += viscosityCoefficient * neighbor->getMass() *
                                      (velocityDiff / neighbor->getDensity()) *
                                      kernelViscosityLaplacian(dist);
                }
            }

            sf::Vector2f gravityForce = gravity * particle->getMass();
            sf::Vector2f acceleration = (pressureForce + viscosityForce + gravityForce) / particle->getMass();
            particle->setAcceleration(acceleration);
        }
    }

    void SPHPhysics::integrate(const std::vector<Particle *> &particles, float dt)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < particles.size(); ++i)
        {
            auto *particle = particles[i];
            sf::Vector2f velocity = particle->getVelocity() + particle->getAcceleration() * dt;
            particle->setVelocity(velocity);
            particle->setPosition(particle->getPosition() + velocity * dt);
        }
    }

    void SPHPhysics::resolveCollisions(const std::vector<Particle*>& particles, Grid* grid, float width, float height) {
		constexpr float PARTICLE_RADIUS = 4.0f;
		constexpr float COLLISION_DAMPING = 0.5f;
		const float minDist = PARTICLE_RADIUS * 2.0f;
		const float minDistSq = minDist * minDist;
		float restitution = 0.5f;
	
		// Compute overlap severity and density metrics
		float maxOverlap = 0.0f;
		float avgNeighbors = 0.0f;
		int maxNeighbors = 0;
		for (auto* p : particles) {
			int neighborCount = p->cachedNeighbors.size();
			avgNeighbors += neighborCount;
			maxNeighbors = std::max(maxNeighbors, neighborCount);
			for (auto* p2 : p->cachedNeighbors) {
				if (p == p2) continue;
				sf::Vector2f delta = p->getPosition() - p2->getPosition();
				float distSq = delta.x * delta.x + delta.y * delta.y;
				if (distSq < minDistSq && distSq > 1e-6f) {
					float dist = std::sqrt(distSq);
					float overlap = minDist - dist;
					maxOverlap = std::max(maxOverlap, overlap);
				}
			}
		}
		avgNeighbors /= particles.size();
	
		// Adaptive maxIterations: base + density + overlap severity
		int maxIterations = 3 + static_cast<int>(std::ceil(avgNeighbors / 4.0f) + maxOverlap / PARTICLE_RADIUS);
		maxIterations = std::min(maxIterations, 15); // Higher cap for dense regions
		std::cout << "Adaptive max iterations: " << maxIterations 
				  << " (avg neighbors: " << avgNeighbors 
				  << ", max neighbors: " << maxNeighbors 
				  << ", max overlap: " << maxOverlap << ")" << std::endl;
	
		// Boundary collisions
		#pragma omp parallel for
		for (size_t i = 0; i < particles.size(); ++i) {
			auto* particle = particles[i];
			sf::Vector2f pos = particle->getPosition();
			sf::Vector2f vel = particle->getVelocity();
			if (pos.x < PARTICLE_RADIUS) { pos.x = PARTICLE_RADIUS; vel.x = -vel.x * boundaryDamping; }
			if (pos.x > width - PARTICLE_RADIUS) { pos.x = width - PARTICLE_RADIUS; vel.x = -vel.x * boundaryDamping; }
			if (pos.y < PARTICLE_RADIUS) { pos.y = PARTICLE_RADIUS; vel.y = -vel.y * boundaryDamping; }
			if (pos.y > height - PARTICLE_RADIUS) { pos.y = height - PARTICLE_RADIUS; vel.y = -vel.y * boundaryDamping; }
			particle->setPosition(pos);
			particle->setVelocity(vel);
		}
	
		// Particle-particle collisions with adaptive iterations
		for (int iter = 0; iter < maxIterations; ++iter) {
			std::vector<std::tuple<Particle*, sf::Vector2f, sf::Vector2f>> globalUpdates;
			#pragma omp parallel 
			{
				std::vector<std::tuple<Particle*, sf::Vector2f, sf::Vector2f>> localUpdates;
				#pragma omp for nowait
				for (size_t i = 0; i < particles.size(); ++i) {
					auto* p1 = particles[i];
					sf::Vector2f pos1 = p1->getPosition();
					for (auto* p2 : p1->cachedNeighbors) {
						if (p1 == p2) continue;
						sf::Vector2f pos2 = p2->getPosition();
						sf::Vector2f delta = pos1 - pos2;
						float distSq = delta.x * delta.x + delta.y * delta.y;
						if (distSq < minDistSq && distSq > 1e-6f) {
							float dist = std::sqrt(distSq);
							float penetration = minDist - dist;
							sf::Vector2f normal = delta / dist;
							// Stronger correction inspired by Simulator
							sf::Vector2f correction = normal * (penetration * 0.25f); // Match Simulator's 0.25f delta
							sf::Vector2f v1 = p1->getVelocity();
							sf::Vector2f v2 = p2->getVelocity();
							sf::Vector2f relVel = v1 - v2;
							float velAlongNormal = relVel.x * normal.x + relVel.y * normal.y;
							sf::Vector2f impulse = {0.0f, 0.0f};
							if (velAlongNormal < 0) {
								float impulseMagnitude = -(1.0f) * velAlongNormal / 2.0f;
								impulse = normal * impulseMagnitude;
							}
							impulse += normal * (penetration * 0.05f); // Fixed correction factor
							localUpdates.push_back(std::make_tuple(p1, pos1 + correction, v1 + impulse));
							localUpdates.push_back(std::make_tuple(p2, pos2 - correction, v2 - impulse));
						}
					}
				}
				#pragma omp critical
				globalUpdates.insert(globalUpdates.end(), localUpdates.begin(), localUpdates.end());
			}
			for (const auto& update : globalUpdates) {
				Particle* p = std::get<0>(update);
				p->setPosition(std::get<1>(update));
				p->setVelocity(std::get<2>(update));
			}
		}
	}

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
}