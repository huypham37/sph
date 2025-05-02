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
    SPHPhysics::SPHPhysics(){
        h = Config::SMOOTHING_RADIUS;
        h2 = h * h;
        viscosityCoefficient = Config::VISCOSITY;
        gasConstant = Config::GAS_CONSTANT;
        restDensity = Config::REST_DENSITY;
        boundaryDamping = Config::BOUNDARY_DAMPING;
        gamma = Config::GAMMA;
        timeStepCounter = 0;
    };

    void SPHPhysics::setSmoothingRadius(float smoothingRadius)
    {
        h = smoothingRadius;
        h2 = h * h;
    }

    void SPHPhysics::computeDensityPressure(const std::vector<Particle *> &particles, Grid *grid)
    {
        // Increment time step counter at the beginning of each simulation step
        timeStepCounter++;
        bool shouldPrintDebug = (timeStepCounter % 1000 == 0);

        // Pre-compute kernel coefficient for optimization
        const float POLY6_COEFF = 4.0f / (M_PI * std::pow(h, 8));

        // Debug - select particles near the bottom of the screen
        debugParticles.clear();
        if (shouldPrintDebug && !particles.empty())
        {
            // Filter particles with y position greater than 490
            for (size_t i = 0; i < particles.size(); ++i)
            {
                auto *particle = particles[i];
                sf::Vector2f position = particle->getPosition();

                if (position.y > 490.0f)
                {
                    debugParticles.insert(i);
                }
            }

            // Output how many particles we're debugging
            if (!debugParticles.empty())
            {
                std::cout << "Debugging " << debugParticles.size() << " particles near bottom" << std::endl;
            }
        }

        for (size_t i = 0; i < particles.size(); ++i)
        {
            auto *particle = particles[i];
            float density = 0.0f;

            // Include self-contribution for density
            density += particle->getMass() * POLY6_COEFF * std::pow(h2, 3);
            // Get neighbors for efficiency in later calculations

            // Sum contribution from neighbors
            for (auto *neighbor : particle->cachedNeighbors)
            {
                if (particle == neighbor)
                    continue;

                sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
                float distSqr = r.x * r.x + r.y * r.y;

                if (distSqr < h2)
                {
                    float h2_r2 = h2 - distSqr;
                    density += neighbor->getMass() * POLY6_COEFF * h2_r2 * h2_r2 * h2_r2;
                }
            }

            // Ensure minimum density
            density = std::max(density, 0.9f * restDensity);
            particle->setDensity(density);

            float density_ratio = density / restDensity;
            // Basic EOS (Equation of State)
            // float pressure = gasConstant * (density - restDensity);
            float pressure_term = std::pow(density_ratio, gamma) - 1.0f;
            float pressure = (this->gasConstant * restDensity / gamma) * pressure_term;
            pressure = std::max(0.0f, pressure);
            particle->setPressure(pressure);

            // Update the debug output code to respect the counter
            if (shouldPrintDebug && debugParticles.find(i) != debugParticles.end())
            {
                int realNeighborCount = 0;
                for (auto *neighbor : particle->cachedNeighbors)
                {
                    if (particle == neighbor)
                        continue;
                    sf::Vector2f r = particle->getPosition() - neighbor->getPosition();
                    float distSqr = r.x * r.x + r.y * r.y;
                    if (distSqr < 4 * h2)
                        realNeighborCount++;
                }

                sf::Vector2f pos = particle->getPosition();
                sf::Vector2f vel = particle->getVelocity();
                std::cout << "Step " << timeStepCounter << " - Particle " << i << " at (" << pos.x << ", " << pos.y << "): "
                          << "Mass: " << particle->getMass() << ", "
                          << "Grid neighbors: " << particle->cachedNeighbors.size()
                          << ", h-radius neighbors: " << realNeighborCount
                          << ", Density: " << particle->getDensity()
                          << ", Pressure: " << particle->getPressure()
                          << ", Velocity: (" << vel.x << ", " << vel.y << ")" << std::endl;
            }
        }
    }

    // Assume Constructor, computeDensityPressure, integrate, resolveCollisions etc. exist

    void SPHPhysics::computeForces(const std::vector<Particle *> &particles, Grid *grid)
    {
        // Kernel gradient constant (for Pressure AND Artificial Viscosity)
        // Note: Using 45 / (pi * h^6) for 3D Spiky gradient, adjust if needed for 2D
        // A common 2D Spiky gradient constant is -30 / (pi * h^5)
        // Let's stick with your value for now, assuming it's correct for your 2D kernel
        const float SPIKY_GRAD_COEFF = -10.0f / (M_PI * std::pow(h, 5));

        // Small value to prevent division by zero and instability
        const float EPS = 1e-6f;
        const float h2 = h * h; // Ensure h2 is accessible (member or calculated)

        // --- Artificial Viscosity Parameters ---
        const float alpha = 1.0f; // TUNABLE: Start 0.1, increase (0.3, 0.5, 1.0) if unstable/energetic
        // const float beta = 0.0f;  // Often zero for liquids
        const float artificial_epsilon = 0.01f * h2; // Small term added to denominator: epsilon*h^2

        // Debugging flag (assuming timeStepCounter is a member variable)
        // bool shouldPrintDebug = (timeStepCounter % 1000 == 0);

        // #pragma omp parallel for // Can parallelize the outer loop
        for (size_t i = 0; i < particles.size(); ++i)
        {
            auto *particle = particles[i];
            sf::Vector2f pressureAcceleration(0.0f, 0.0f);
            sf::Vector2f viscosityAcceleration(0.0f, 0.0f); // Accumulates artificial viscosity acceleration

            // Get particle properties once
            sf::Vector2f pos_i = particle->getPosition();
            sf::Vector2f vel_i = particle->getVelocity();
            float density_i = particle->getDensity();
            float pressure_i = particle->getPressure();

            // Make sure density isn't too close to zero before division
            if (density_i < EPS)
                density_i = EPS; // Safeguard

            for (auto *neighbor : particle->cachedNeighbors)
            {
                if (particle == neighbor)
                    continue;

                // Get neighbor properties once
                sf::Vector2f pos_j = neighbor->getPosition();
                sf::Vector2f vel_j = neighbor->getVelocity();
                float density_j = neighbor->getDensity();
                float pressure_j = neighbor->getPressure();
                float mass_j = neighbor->getMass();

                // Make sure density isn't too close to zero before division
                if (density_j < EPS)
                    density_j = EPS; // Safeguard

                sf::Vector2f r_ij = pos_i - pos_j; // Vector from j to i
                float distSqr = r_ij.x * r_ij.x + r_ij.y * r_ij.y;

                if (distSqr < h2 && distSqr > EPS) // Check within h and avoid self/coincident
                {
                    float dist = std::sqrt(distSqr);
                    sf::Vector2f dir_ij = r_ij / dist; // Normalized direction from j to i

                    // --- Pressure Acceleration Calculation (Unchanged) ---
                    float pressureTerm = -(pressure_i / (density_i * density_i) +
                                           pressure_j / (density_j * density_j));

                    float h_r = h - dist;
                    // Use Spiky gradient for pressure too (common practice)
                    sf::Vector2f gradW_spiky = SPIKY_GRAD_COEFF * h_r * h_r * dir_ij;
                    pressureAcceleration += mass_j * pressureTerm * gradW_spiky;

                    // --- Monaghan Artificial Viscosity Acceleration Calculation ---
                    sf::Vector2f vel_ij = vel_i - vel_j;                   // Relative velocity v_i - v_j
                    float v_dot_r = vel_ij.x * r_ij.x + vel_ij.y * r_ij.y; // v_ij dot r_ij

                    // Only apply viscosity if particles are *approaching* each other
                    if (v_dot_r < 0.0f)
                    {
                        // Average density
                        float rho_avg = 0.5f * (density_i + density_j);

                        // Speed of sound estimate (can be simplified or made more complex)
                        // Using sqrt(kappa/rho) is more related to Tait EOS.
                        // A simpler estimate often used here:
                        float c_i = std::sqrt(std::max(0.0f, pressure_i / density_i)); // Speed of sound estimate at i
                        float c_j = std::sqrt(std::max(0.0f, pressure_j / density_j)); // Speed of sound estimate at j
                        float c_avg = 0.5f * (c_i + c_j);

                        // Mu term from Monaghan's paper
                        float mu_ij = (h * v_dot_r) / (distSqr + artificial_epsilon); // distSqr = r_ij^2

                        // Viscosity term (PI_ij) - Ignoring beta term for liquids
                        float PI_ij = (-alpha * c_avg * mu_ij) / rho_avg;

                        // Viscosity Acceleration contribution = - sum [ m_j * PI_ij * gradW_ij ]
                        // Note: gradW_spiky is already computed above for pressure
                        viscosityAcceleration += -mass_j * PI_ij * gradW_spiky;
                    }
                    // --- End Artificial Viscosity ---
                }
            } // End neighbor loop

            // Set total acceleration = Pressure Acceleration + Viscosity Acceleration + Gravity
            particle->setAcceleration(pressureAcceleration + viscosityAcceleration + gravity);
        } // End particle loop
    } // End computeForces

    void SPHPhysics::integrate(const std::vector<Particle *> &particles, float dt)
    {
        // #pragma omp parallel for
        for (size_t i = 0; i < particles.size(); ++i)
        {
            auto *particle = particles[i];
            sf::Vector2f velocity = particle->getVelocity() + particle->getAcceleration() * dt;
            particle->setVelocity(velocity);
            particle->setPosition(particle->getPosition() + velocity * dt);
        }
    }

    // void SPHPhysics::integrateVerlet(const std::vector<Particle *> &particles, float dt, Grid* grid, float width, float height) {
    //     // First half of Verlet using existing forces
    //     for (auto* particle: particles)
    //     {
    //         sf::Vector2f halfVelocity = particle->getVelocity() +
    //                                     0.5f * dt * particle->getAcceleration();
    //         // update position
    //         sf::Vector2f newPosition = particle->getPosition() + halfVelocity * dt;
    //         particle->setPosition(newPosition);
    //         particle->setImmediateVelocity(halfVelocity);
    //     }

    //     // Update grid with new positions:
    //     particles->updateGrid(particles);

    //     computeDensityPressure(particles, grid);

    //     computeForces(particles, grid);

    //     for (auto* particle: particles)
    //     {
    //         sf::Vector2f finalVelocity = particle->getImmediateVelocity() +
    //                                      particle->getAcceleration() * (0.5f * dt);
    //         particle->setVelocity(finalVelocity);
    //     }

    //     resolveCollisions(particles, grid, width, height);
    // }

    void SPHPhysics::computeBoundaryForces(const std::vector<Particle *> &particles, float width, float height)
    {
        constexpr float PARTICLE_RADIUS = 5.0f;

        // Boundary parameters - adjusted for mass = 100.0f
        const float boundaryDistance = 1.5f * PARTICLE_RADIUS; // Detection distance
        const float boundaryStiffness = 10000.0f;             // Adjusted for mass=100.0f
        const float boundaryDecay = 2.0f;                      // How quickly force decays with distance

#pragma omp parallel for
        for (size_t i = 0; i < particles.size(); ++i)
        {
            auto *particle = particles[i];
            sf::Vector2f pos = particle->getPosition();
            sf::Vector2f boundaryForce(0.0f, 0.0f);

            // Left wall repulsion
            if (pos.x < boundaryDistance)
            {
                float dist = pos.x;
                float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
                float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
                boundaryForce.x += forceMagnitude;
            }

            // Right wall repulsion
            if (pos.x > width - boundaryDistance)
            {
                float dist = width - pos.x;
                float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
                float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
                boundaryForce.x -= forceMagnitude;
            }

            // Top wall repulsion
            if (pos.y < boundaryDistance)
            {
                float dist = pos.y;
                float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
                float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
                boundaryForce.y += forceMagnitude;
            }

            // Bottom wall repulsion
            if (pos.y > height - boundaryDistance)
            {
                float dist = height - pos.y;
                float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
                float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
                boundaryForce.y -= forceMagnitude;
            }

            // Add boundary acceleration to existing acceleration
            sf::Vector2f currentAcc = particle->getAcceleration();
            sf::Vector2f boundaryAcc = boundaryForce / particle->getDensity();
            particle->setAcceleration(currentAcc + boundaryAcc);
        }
    }

    void SPHPhysics::resolveCollisions(const std::vector<Particle *> &particles, Grid *grid, float width, float height)
    {
        constexpr float PARTICLE_RADIUS = 5.0f;

#pragma omp parallel for
        for (size_t i = 0; i < particles.size(); ++i)
        {
            auto *particle = particles[i];
            sf::Vector2f pos = particle->getPosition();
            sf::Vector2f vel = particle->getVelocity();

            // Simple boundary conditions with damping - failsafe only
            if (pos.x < PARTICLE_RADIUS)
            {
                pos.x = PARTICLE_RADIUS;
                vel.x = -vel.x * boundaryDamping;
            }
            if (pos.x > width - PARTICLE_RADIUS)
            {
                pos.x = width - PARTICLE_RADIUS;
                vel.x = -vel.x * boundaryDamping;
            }
            if (pos.y < PARTICLE_RADIUS)
            {
                pos.y = PARTICLE_RADIUS;
                vel.y = -vel.y * boundaryDamping;
            }
            if (pos.y > height - PARTICLE_RADIUS)
            {
                pos.y = height - PARTICLE_RADIUS;
                vel.y = -vel.y * boundaryDamping;
            }

            particle->setPosition(pos);
            particle->setVelocity(vel);
        }
    }
}
