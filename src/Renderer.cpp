#include "Renderer.hpp"
#include <iostream>
#include <algorithm>

namespace sph
{

	Renderer::Renderer()
		: fontLoaded(false)
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

} // namespace sph