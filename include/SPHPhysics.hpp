#pragma once

#include "Grid.hpp" // Added missing include for Grid class
#include "Particle.hpp"
#include "SPHConfig.hpp"

#include <SFML/System/Vector2.hpp>
#include <unordered_set> // Added missing include for unordered_set
#include <vector>

namespace sph {

/**
 * @brief Handles SPH physics calculations
 *
 * Responsible for implementing the SPH fluid simulation algorithm,
 * including density, pressure, forces, and physical interactions.
 */
class SPHPhysics {
public:
  /**
   * @brief Constructor
   */
  SPHPhysics();
  /**
   * @brief Set simulation parameters
   *
   * @param smoothingRadius SPH kernel smoothing radius
   */
  void setSmoothingRadius(float smoothingRadius);

  /**
   * @brief Set gravity force
   *
   * @param x X component of gravity
   * @param y Y component of gravity
   */
  void setGravity(float x, float y) { gravity = {x, y}; }

  /**
   * @brief Set viscosity coefficient
   *
   * @param viscosity Viscosity coefficient value
   */
  void setViscosity(float viscosity) { viscosityCoefficient = viscosity; }

  /**
   * @brief Set gas constant for pressure calculation
   *
   * @param gasConstantValue Gas constant value
   */
  void setGasConstant(float gasConstantValue) {
    gasConstant = gasConstantValue;
  }

  /**
   * @brief Set rest density for pressure calculation
   *
   * @param density Rest density value
   */
  void setRestDensity(float density) { restDensity = density; }

  /**
   * @brief Set boundary damping factor
   *
   * @param damping Boundary damping coefficient
   */
  void setBoundaryDamping(float damping) { boundaryDamping = damping; }

  /**
   * @brief Compute density and pressure for each particle
   *
   * @param particles Vector of particles to process
   * @param grid Spatial partitioning grid for neighbor search
   */
  void computeDensityPressure(const std::vector<Particle *> &particles,
                              Grid *grid);

  /**
   * @brief Compute forces (pressure, viscosity, gravity) for each particle
   *
   * @param particles Vector of particles to process
   * @param grid Spatial partitioning grid for neighbor search
   */
  void computeForces(const std::vector<Particle *> &particles, Grid *grid);

  /**
   * @brief Integrate particle positions and velocities
   *
   * @param particles Vector of particles to process
   * @param dt Time step for integration
   */
  void integrate(const std::vector<Particle *> &particles, float dt);

  /**
   * @brief Resolve collisions between particles and boundaries
   *
   * @param particles Vector of particles to process
   * @param grid Spatial partitioning grid for neighbor search
   * @param width Width of simulation domain
   * @param height Height of simulation domain
   */
  void resolveCollisions(const std::vector<Particle *> &particles, Grid *grid,
                         float width, float height);

  /**
   * @brief Compute forces exerted by boundaries on particles
   *
   * @param particles Vector of particles to process
   * @param width Width of simulation domain
   * @param height Height of simulation domain
   */
  void computeBoundaryForces(const std::vector<Particle *> &particles,
                             float width, float height);

  /**
   * @brief Compute density and pressure for a single particle
   *
   * @param particle Pointer to the particle to process
   */
  void computeDensityPressureForParticle(Particle* particle);

  /**
   * @brief Compute forces (pressure, viscosity, gravity) for a single particle
   *
   * @param particle Pointer to the particle to process
   */
  void computeForcesForParticle(Particle* particle);

  /**
   * @brief Compute boundary forces for a single particle
   *
   * @param particle Pointer to the particle to process
   * @param width Width of simulation domain
   * @param height Height of simulation domain
   */
  void computeBoundaryForcesForParticle(Particle* particle, float width, float height);

  /**
   * @brief Integrate position and velocity for a single particle
   *
   * @param particle Pointer to the particle to process
   * @param dt Time step for integration
   */
  void integrateParticle(Particle* particle, float dt);

  /**
   * @brief Resolve collisions for a single particle
   *
   * @param particle Pointer to the particle to process
   * @param width Width of simulation domain
   * @param height Height of simulation domain
   */
  void resolveCollisionsForParticle(Particle* particle, float width, float height);

  // Getter methods for SPH parameters
  float getSmoothingRadius() const { return h; }
  float getViscosity() const { return viscosityCoefficient; }
  float getGasConstant() const { return gasConstant; }
  float getRestDensity() const { return restDensity; }
  sf::Vector2f getGravity() const { return gravity; }

  // Expose kernel functions for parallel tasks
  float kernelPoly6(float distSquared);
  sf::Vector2f kernelGradSpiky(float distance, const sf::Vector2f &direction);
  float kernelViscosityLaplacian(float distance);

private:
  // SPH Parameters
  sf::Vector2f gravity;                      // Gravity force
  float h;                                   // Smoothing radius
  float h2;                                  // h^2 (pre-computed)
  float viscosityCoefficient;                // Viscosity coefficient
  float gasConstant;                         // Gas constant for pressure
  float restDensity;                         // Rest density
  float boundaryDamping;                     // Boundary collision damping
  float gamma;                               // For equation of state
  int timeStepCounter;                       // For debugging
  std::unordered_set<size_t> debugParticles; // For tracking specific particles
};

} // namespace sph
