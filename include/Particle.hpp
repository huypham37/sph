#pragma once

#include <SFML/Graphics.hpp>

class Particle
{
public:
	Particle(float x, float y);
	~Particle() = default;
	
	std::vector<Particle*> cachedNeighbors;
	void update(float dt);
	// Update the particle color based on velocity
	void updateColor();
	void updateVisuals();
	void draw(sf::RenderWindow &window);

	sf::Vector2f getPosition() const { return position; }
	sf::Vector2f getVelocity() const { return velocity; }
	sf::Vector2f getAcceleration() const { return acceleration; }
	float getMass() const { return mass; }
	float getDensity() const { return density; }
	float getPressure() const { return pressure; }

	void setPosition(const sf::Vector2f &pos) { position = pos; }
	void setVelocity(const sf::Vector2f &vel) { velocity = vel; }
	void setAcceleration(const sf::Vector2f &acc) { acceleration = acc; }
	void setDensity(float d) { density = d; }
	void setPressure(float p) { pressure = p; }

private:
	sf::Vector2f position;
	sf::Vector2f velocity;
	sf::Vector2f acceleration;
	float mass;
	float density;
	float pressure;
	sf::CircleShape shape;
	sf::Color baseColor; // Base color for the particle

	static constexpr float RADIUS = 4.0f;
	static constexpr float DEFAULT_MASS = 1.0f;
	
};