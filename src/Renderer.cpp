#include "Renderer.hpp"
#include <iostream>
#include <algorithm>
#include <numeric>

namespace sph
{

	Renderer::Renderer()
		: fontLoaded(false),
		  visualizeSubdomains(false),
		  visualizeLoadBalance(false)
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

	void Renderer::drawSubdomainLoadInfo(
		const std::vector<std::unique_ptr<parallel::Subdomain>> &subdomains,
		sf::RenderWindow &window)
	{
		if (!visualizeLoadBalance || subdomains.empty() || !fontLoaded)
		{
			return;
		}

		// Calculate statistics for normalization
		float totalTime = 0.0f;
		float maxTime = 0.0f;
		std::vector<float> times;

		for (const auto &subdomain : subdomains)
		{
			float time = subdomain->getLastComputationTime();
			totalTime += time;
			maxTime = std::max(maxTime, time);
			times.push_back(time);
		}

		float avgTime = totalTime / subdomains.size();

		// Draw load info for each subdomain
		const float barHeight = 10.0f;
		const float barWidth = 60.0f;
		const float padding = 5.0f;

		// Position variables - top right corner
		float xPos = window.getSize().x - 270.0f; // Position from right side
		float yPos = 10.0f;						  // Starting from top

		// Draw title
		sf::Text titleText(font, "Load Balancing Metrics", 14);
		titleText.setFillColor(sf::Color::White);
		titleText.setPosition({xPos, yPos});
		window.draw(titleText);
		yPos += 25.0f;

		// Draw average computation time info
		sf::Text avgText(font, "Avg time: " + std::to_string(avgTime * 1000.0f) + " ms", 12);
		avgText.setFillColor(sf::Color::White);
		avgText.setPosition({xPos, yPos});
		window.draw(avgText);
		yPos += 20.0f;

		// Draw load bars for each subdomain
		for (size_t i = 0; i < subdomains.size(); i++)
		{
			const auto &subdomain = subdomains[i];
			float time = subdomain->getLastComputationTime();
			float normalizedTime = (maxTime > 0.0f) ? (time / maxTime) : 0.0f;

			// Draw subdomain ID and stats
			std::string info = "Domain " + std::to_string(subdomain->getId()) +
							   ": " + std::to_string(subdomain->getParticles().size()) + " particles, " +
							   std::to_string(time * 1000.0f) + " ms";

			sf::Text text(font, info, 12);
			text.setFillColor(sf::Color::White);
			text.setPosition({xPos, yPos});
			window.draw(text);

			// Draw load bar background
			sf::RectangleShape barBg;
			barBg.setPosition({xPos + 10.0f, yPos + padding + 15.0f});
			barBg.setSize({barWidth, barHeight});
			barBg.setFillColor(sf::Color(50, 50, 50));
			window.draw(barBg);

			// Draw actual load bar
			sf::RectangleShape bar;
			bar.setPosition({xPos + 10.0f, yPos + padding + 15.0f});
			bar.setSize({barWidth * normalizedTime, barHeight});
			bar.setFillColor(getLoadColor(normalizedTime));
			window.draw(bar);

			// Color the subdomain based on load
			if (visualizeSubdomains)
			{
				sf::RectangleShape rect;
				rect.setPosition({subdomain->getX(), subdomain->getY()});
				rect.setSize({subdomain->getWidth(), subdomain->getHeight()});
				rect.setFillColor(sf::Color(getLoadColor(normalizedTime).r,
											getLoadColor(normalizedTime).g,
											getLoadColor(normalizedTime).b, 40));
				rect.setOutlineColor(sf::Color(255, 255, 255, 80));
				rect.setOutlineThickness(1.0f);
				window.draw(rect);
			}

			yPos += barHeight + padding + 30.0f;
		}
	}

	sf::Color Renderer::getLoadColor(float load) const
	{
		// Green (low load) to Yellow to Red (high load)
		if (load < 0.5f)
		{
			// Green to Yellow
			int r = static_cast<int>(255 * (load * 2.0f));
			return sf::Color(r, 255, 0);
		}
		else
		{
			// Yellow to Red
			int g = static_cast<int>(255 * (2.0f - load * 2.0f));
			return sf::Color(255, g, 0);
		}
	}

} // namespace sph