#include "Particle.hpp"
#include <cmath>

Particle::Particle(float x, float y)
	: position({x, y}),
	  velocity({0.0f, 0.0f}),
	  acceleration({0.0f, 0.0f}),
	  mass(DEFAULT_MASS),
	  density(0.0f),
	  pressure(0.0f),
	  shape(RADIUS)
{
	// Base color for still particles
	baseColor = sf::Color(40, 120, 255, 220);
	shape.setFillColor(baseColor);
	shape.setOrigin({RADIUS, RADIUS});
	shape.setPosition(position);

	// Add slight outline for better visibility
	shape.setOutlineThickness(0.5f);
	shape.setOutlineColor(sf::Color(10, 50, 150, 100));
}

void Particle::update(float dt)
{
	// Update position
	position += velocity * dt;

	// Update color based on velocity
	updateColor();

	// Update shape position for rendering
	shape.setPosition(position);
}

void Particle::draw(sf::RenderWindow &window)
{
	window.draw(shape);
}

void Particle::updateColor()
{
	// Calculate speed of the particle
	float speed = std::sqrt(velocity.x * velocity.x + velocity.y * velocity.y);

	// Map speed to a color gradient
	// - Low speed: blue (cold)
	// - Medium speed: light blue/cyan
	// - High speed: white with blue tint

	constexpr float MAX_SPEED = 100.0f;
	float speedFactor = std::min(speed / MAX_SPEED, 1.0f);

	// Blue to white gradient
	int r = baseColor.r + static_cast<int>((255 - baseColor.r) * speedFactor * 0.8f);
	int g = baseColor.g + static_cast<int>((255 - baseColor.g) * speedFactor);
	int b = baseColor.b + static_cast<int>((255 - baseColor.b) * speedFactor * 0.4f);

	// Slight glow effect for fast-moving particles
	int alpha = baseColor.a + static_cast<int>((255 - baseColor.a) * speedFactor * 0.3f);

	shape.setFillColor(sf::Color(r, g, b, alpha));

	// Adjust radius slightly based on speed for more dynamic appearance
	float radiusMultiplier = 1.0f + speedFactor * 0.15f;
	shape.setRadius(RADIUS * radiusMultiplier);
	shape.setOrigin({RADIUS * radiusMultiplier, RADIUS * radiusMultiplier});

	// Update outline based on speed as well
	if (speedFactor > 0.7f)
	{
		shape.setOutlineColor(sf::Color(200, 230, 255, 120));
		shape.setOutlineThickness(0.8f);
	}
	else
	{
		shape.setOutlineColor(sf::Color(10, 50, 150, 100));
		shape.setOutlineThickness(0.5f);
	}
}