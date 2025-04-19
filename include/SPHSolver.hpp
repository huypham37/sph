#pragma once

#include <vector>
#include <memory>
#include "Particle.hpp"
#include "Grid.hpp"
#include "SPHPhysics.hpp" // Include SPHPhysics header

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

	// Dynamic particle management
	void addParticles(int count);
	void removeParticles(int count);

	// SPH parameters setters (delegate to physics engine)
	void setGravity(float x, float y) { physicsEngine.setGravity(x, y); } // Fixed signature
	void setViscosity(float v) { physicsEngine.setViscosity(v); }
	void setGasConstant(float k) { physicsEngine.setGasConstant(k); }
	void setRestDensity(float d) { physicsEngine.setRestDensity(d); }
	void setSmoothingRadius(float radius);

private:
	std::vector<Particle *> particles;
	Grid *grid;
	sph::SPHPhysics physicsEngine;

	float width;
	float height;

	// Font for visualization text (if needed)
	sf::Font font;
	bool fontLoaded;
};