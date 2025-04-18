#pragma once

#include "SPHSolver.hpp"
#include "parallel/DomainDecomposer.hpp"
#include "parallel/BoundaryManager.hpp"
#include "parallel/LoadBalancer.hpp"
#include <memory>

namespace sph
{
	namespace parallel
	{

		/**
		 * @brief Parallel implementation of the SPH solver using domain decomposition
		 *
		 * This class extends the basic SPHSolver to use domain decomposition for
		 * parallel processing across multiple threads.
		 */
		class ParallelSPHSolver : public SPHSolver
		{
		public:
			/**
			 * @brief Constructor
			 *
			 * @param width Width of the simulation domain
			 * @param height Height of the simulation domain
			 * @param numThreads Number of threads to use (defaults to max available)
			 */
			ParallelSPHSolver(float width, float height, int numThreads = 0);

			/**
			 * @brief Destructor
			 */
			~ParallelSPHSolver() override;

			/**
			 * @brief Update the simulation state
			 *
			 * @param dt Time step
			 */
			void update(float dt) override;

			/**
			 * @brief Draw the simulation
			 *
			 * @param window SFML window to draw to
			 */
			void draw(sf::RenderWindow &window) override;

			/**
			 * @brief Set the domain decomposer to use
			 *
			 * @param decomposer Domain decomposition strategy
			 */
			void setDomainDecomposer(std::unique_ptr<DomainDecomposer> decomposer);

			/**
			 * @brief Set the load balancer to use
			 *
			 * @param balancer Load balancing strategy
			 */
			void setLoadBalancer(std::unique_ptr<LoadBalancer> balancer);

			/**
			 * @brief Get the number of subdomains (approximately equals thread count)
			 *
			 * @return Current number of subdomains
			 */
			int getNumSubdomains() const;

			/**
			 * @brief Set the number of subdomains/threads to use
			 *
			 * @param count Number of subdomains (0 for auto)
			 */
			void setNumSubdomains(int count);

			/**
			 * @brief Toggle visualization of subdomains
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
			// Domain decomposition components
			std::unique_ptr<DomainDecomposer> domainDecomposer;
			std::unique_ptr<BoundaryManager> boundaryManager;
			std::unique_ptr<LoadBalancer> loadBalancer;

			// Subdomains created by the decomposer
			std::vector<std::unique_ptr<Subdomain>> subdomains;

			// Thread management
			int numSubdomains;

			// Visualization options
			bool visualizeSubdomains;

			/**
			 * @brief Initialize domain decomposition
			 */
			void initializeDecomposition();

			/**
			 * @brief Compute density and pressure in parallel across subdomains
			 */
			void computeDensityPressureParallel();

			/**
			 * @brief Compute forces in parallel across subdomains
			 */
			void computeForcesParallel();

			/**
			 * @brief Integrate particles in parallel across subdomains
			 *
			 * @param dt Time step
			 */
			void integrateParallel(float dt);

			/**
			 * @brief Update load balancing based on subdomain performance
			 */
			void updateLoadBalancing();

			/**
			 * @brief Draw subdomain boundaries for visualization
			 *
			 * @param window SFML window to draw to
			 */
			void drawSubdomains(sf::RenderWindow &window);
		};

	} // namespace parallel
} // namespace sph