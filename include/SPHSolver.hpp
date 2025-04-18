#pragma once

#include <vector>
#include <memory>
#include "Particle.hpp"
#include "Grid.hpp"
// Include parallel components
#include "parallel/Subdomain.hpp"
#include "parallel/GridDomainDecomposer.hpp"
#include "parallel/BoundaryManager.hpp"

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

	// Parallel options
	void setParallelizationEnabled(bool enabled);
	bool isParallelizationEnabled() const { return parallelizationEnabled; }
	void setVisualizeSubdomains(bool enabled) { visualizeSubdomains = enabled; }
	bool isVisualizeSubdomains() const { return visualizeSubdomains; }

private:
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

	// Parallel components
	bool parallelizationEnabled;
	bool visualizeSubdomains;
	std::unique_ptr<sph::parallel::GridDomainDecomposer> domainDecomposer;
	std::unique_ptr<sph::parallel::BoundaryManager> boundaryManager;
	std::vector<std::unique_ptr<sph::parallel::Subdomain>> subdomains;

	// Font for visualization text
	sf::Font font;
	bool fontLoaded;

	// Initialize domain decomposition
	void initializeParallelComponents();
	void updateDomainDecomposition();

	// SPH Methods
	void computeDensityPressure();
	void computeForces();
	void integrate(float dt);
	void resolveCollisions();

	// Parallel SPH methods
	void computeDensityPressureParallel();
	void computeForcesParallel();
	void integrateParallel(float dt);
	void drawSubdomains(sf::RenderWindow &window);

	// SPH Kernels
	float kernelPoly6(float distSqr);
	sf::Vector2f kernelGradSpiky(float dist, const sf::Vector2f &dir);
	float kernelViscosityLaplacian(float dist);
};