#pragma once

#include <vector>
#include <unordered_map>
#include "Particle.hpp"

namespace sph
{

	class Grid
	{
	public:
		Grid(float width, float height, float cellSize);
		~Grid() = default;

		void clear();
		void insertParticle(sph::Particle *particle);
		std::vector<sph::Particle *> getNeighbors(float x, float y, float radius);
		std::vector<sph::Particle *> getNeighbors(const sph::Particle *particle, float radius);
		void updateGrid(const std::vector<sph::Particle *> &particles);
		// Add a getter for cell size to make its use explicit
		float getCellSize() const { return cellSize; }

	private:
		float cellSize;
		float invCellSize; // Inverse of cell size for optimization
		int gridWidth;
		int gridHeight;

		// Grid cell to particles mapping
		std::unordered_map<int, std::vector<sph::Particle *>> cells;

		// Helper methods
		int getCellIndex(float x, float y) const;
		int getCellIndex(int cellX, int cellY) const;
		std::pair<int, int> positionToCell(float x, float y) const;
	};

} // namespace sph