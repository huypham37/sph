#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include <memory>
#include "Particle.hpp"
#include "parallel/Subdomain.hpp"

namespace sph
{

	/**
	 * @brief Handles rendering of simulation elements
	 *
	 * Responsible for drawing particles and visualization elements
	 * like domain decomposition boundaries.
	 */
	class Renderer
	{
	public:
		/**
		 * @brief Constructor
		 */
		Renderer();

		/**
		 * @brief Draw all particles
		 *
		 * @param particles Vector of particles to draw
		 * @param window SFML window to draw to
		 */
		void drawParticles(const std::vector<Particle *> &particles, sf::RenderWindow &window);

		/**
		 * @brief Draw subdomain boundaries
		 *
		 * @param subdomains Vector of subdomains to visualize
		 * @param window SFML window to draw to
		 */
		void drawSubdomains(
			const std::vector<std::unique_ptr<parallel::Subdomain>> &subdomains,
			sf::RenderWindow &window);

		/**
		 * @brief Enable or disable subdomain visualization
		 *
		 * @param enabled Whether to draw subdomain boundaries
		 */
		void setVisualizeSubdomains(bool enabled) { visualizeSubdomains = enabled; }

		/**
		 * @brief Check if subdomain visualization is enabled
		 *
		 * @return Current visualization state
		 */
		bool isVisualizeSubdomains() const { return visualizeSubdomains; }

	private:
		sf::Font font;			  // Font for text rendering
		bool fontLoaded;		  // Whether font was successfully loaded
		bool visualizeSubdomains; // Whether to draw subdomain boundaries
	};

} // namespace sph