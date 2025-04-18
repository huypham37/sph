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
		 * @brief Draw subdomain load information
		 *
		 * @param subdomains Vector of subdomains with load info
		 * @param window SFML window to draw to
		 */
		void drawSubdomainLoadInfo(
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

		/**
		 * @brief Enable or disable load balance visualization
		 *
		 * @param enabled Whether to draw load balance info
		 */
		void setVisualizeLoadBalance(bool enabled) { visualizeLoadBalance = enabled; }

		/**
		 * @brief Check if load balance visualization is enabled
		 *
		 * @return Current load balance visualization state
		 */
		bool isVisualizeLoadBalance() const { return visualizeLoadBalance; }

	private:
		sf::Font font;			   // Font for text rendering
		bool fontLoaded;		   // Whether font was successfully loaded
		bool visualizeSubdomains;  // Whether to draw subdomain boundaries
		bool visualizeLoadBalance; // Whether to draw load balance info

		/**
		 * @brief Calculate color representing computational load
		 *
		 * @param load Normalized load value (0.0-1.0)
		 * @return Color ranging from green (low load) to red (high load)
		 */
		sf::Color getLoadColor(float load) const;
	};

} // namespace sph