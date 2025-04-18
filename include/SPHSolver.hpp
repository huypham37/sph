#pragma once

#include <vector>
#include "Particle.hpp"
#include "Grid.hpp"

namespace sph
{
	namespace parallel
	{
		class ParallelSPHSolver; // Forward declaration
	}
}

class SPHSolver
{
public:
	SPHSolver(float width, float height);
	virtual ~SPHSolver(); // Make destructor virtual for proper inheritance

	virtual void update(float dt);				 // Make update virtual for override
	virtual void draw(sf::RenderWindow &window); // Make draw virtual for override
	void addParticle(float x, float y);
	void reset();

	// Initialize simulation with default particles
	void initializeDefaultParticles(int count);

	// Mouse interaction with particles
	void applyMouseForce(const sf::Vector2f &mousePos, float strength);

	// Get particle count
	size_t getParticleCount() const { return particles.size(); }

	// Thread management
	void setNumThreads(int num);
	int getNumThreads() const { return numThreads; }
	int getMaxThreads() const;

	// Dynamic particle management
	void addParticles(int count);
	void removeParticles(int count);

	// SPH parameters
	void setGravity(float x, float y) { gravity = sf::Vector2f{x, y}; } // Updated to use braced initialization
	void setViscosity(float v) { viscosityCoefficient = v; }
	void setGasConstant(float k) { gasConstant = k; }
	void setRestDensity(float d) { restDensity = d; }

	// Access to SPH kernel methods for derived classes
	float kernelPoly6(float distSqr) const;
	sf::Vector2f kernelGradSpiky(float dist, const sf::Vector2f &dir) const;
	float kernelViscosityLaplacian(float dist) const;

protected:
	// Friend declaration for ParallelSPHSolver
	friend class sph::parallel::ParallelSPHSolver;

	std::vector<Particle *> particles;
	Grid *grid;

	float width;
	float height;

	// Thread control
	int numThreads;

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
};