#include "Renderer.hpp"
#include <iostream>

namespace sph
{

	Renderer::Renderer()
		: fontLoaded(false),
		  visualizeSubdomains(false)
	{
		// Try to load a system font for visualization
		// This is optional; we'll handle the case when the font can't be loaded
		fontLoaded = font.openFromFile("/System/Library/Fonts/Helvetica.ttc");
		if (!fontLoaded)
		{
			std::cout << "Warning: Failed to load font for visualization" << std::endl;
		}
	}

	void Renderer::drawParticles(const std::vector<Particle *> &particles, sf::RenderWindow &window)
	{
		for (auto *particle : particles)
		{
			particle->draw(window);
		}
	}

	void Renderer::drawSubdomains(
		const std::vector<std::unique_ptr<parallel::Subdomain>> &subdomains,
		sf::RenderWindow &window)
	{
		if (!visualizeSubdomains || subdomains.empty())
		{
			return;
		}

		// Create rectangles to represent each subdomain
		for (const auto &subdomain : subdomains)
		{
			sf::RectangleShape rect;
			rect.setPosition({subdomain->getX(), subdomain->getY()});	   // Use SFML 3.0 Vector2f syntax
			rect.setSize({subdomain->getWidth(), subdomain->getHeight()}); // Use SFML 3.0 Vector2f syntax
			rect.setFillColor(sf::Color::Transparent);
			rect.setOutlineColor(sf::Color(255, 255, 255, 80)); // Semi-transparent white
			rect.setOutlineThickness(1.0f);

			// Draw the rectangle
			window.draw(rect);

			// Draw subdomain ID only if we have a font
			if (fontLoaded)
			{
				// Proper SFML 3.0 Text construction requires a font
				sf::Text text(font, std::to_string(subdomain->getId()), 12); // Use SFML 3.0 constructor
				text.setFillColor(sf::Color(255, 255, 255, 128));			 // Semi-transparent white

				// Set position with Vector2f syntax
				text.setPosition({subdomain->getX() + subdomain->getWidth() / 2.0f - 5.0f,
								  subdomain->getY() + subdomain->getHeight() / 2.0f - 8.0f});

				window.draw(text);
			}
		}
	}

} // namespace sph