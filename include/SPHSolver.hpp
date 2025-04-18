#pragma once

#include <vector>
#include "Particle.hpp"
#include "Grid.hpp"

class SPHSolver
{
public:
	SPHSolver(float width, float height);
	~SPHSolver();

	void update(float dt);
	void draw(sf::RenderWindow &window);
	void addParticle(float x, float y);
	void reset();

	// Initialize simulation with default particles
	void initializeDefaultParticles(int count);

	// Mouse interaction with particles
	void applyMouseForce(const sf::Vector2f &mousePos, float strength);

	// Get particle count
	size_t getParticleCount() const { return particles.size(); }

	// SPH parameters
	void setGravity(float x, float y) { gravity = sf::Vector2f{x, y}; } // Updated to use braced initialization
	void setViscosity(float v) { viscosityCoefficient = v; }
	void setGasConstant(float k) { gasConstant = k; }
	void setRestDensity(float d) { restDensity = d; }

private:
	std::vector<Particle *> particles;
	Grid *grid;

	float width;
	float height;

	// SPH Parameters
	sf::Vector2f gravity;
	float h;  // Smoothing radius
	float h2; // h^2
	float viscosityCoefficient;
	float gasConstant;
	float restDensity;
	float boundaryDamping;

	// SPH Methods
	void computeDensityPressure();
	void computeForces();
	void integrate(float dt);
	void resolveCollisions();

	// SPH Kernels
	float kernelPoly6(float distSqr);
	sf::Vector2f kernelGradSpiky(float dist, const sf::Vector2f &dir);
	float kernelViscosityLaplacian(float dist);
};