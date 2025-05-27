#include "Grid.hpp"
#include <cmath>

namespace sph
{

  Grid::Grid(float width, float height, float cellSize)
      : cellSize(cellSize), invCellSize(1.0f / cellSize),
        gridWidth(static_cast<int>(std::ceil(width / cellSize))),
        gridHeight(static_cast<int>(std::ceil(height / cellSize))) {}

  void Grid::clear() { cells.clear(); }

  void Grid::insertParticle(sph::Particle *particle)
  {
    auto pos = particle->getPosition();
    int index = getCellIndex(pos.x, pos.y);
    cells[index].push_back(particle);
  }

  void Grid::updateGrid(const std::vector<sph::Particle *> &particles)
  {
    this->clear();

    for (auto *particle : particles)
    {
      insertParticle(particle);
    }
  }

  std::vector<sph::Particle *> Grid::getNeighbors(float x, float y, float radius)
  {
    std::vector<sph::Particle *> neighbors;

    // Calculate cell coordinates and radius in cells
    auto [cellX, cellY] = positionToCell(x, y);
    int cellRadius = static_cast<int>(std::ceil(radius * invCellSize));
    float radiusSquared = radius * radius; // For distance check optimization

    // Iterate through potential neighbor cells
    for (int i = -cellRadius; i <= cellRadius; ++i)
    {
      for (int j = -cellRadius; j <= cellRadius; ++j)
      {
        int neighborCellX = cellX + i;
        int neighborCellY = cellY + j;

        // Skip if outside grid bounds
        if (neighborCellX < 0 || neighborCellX >= gridWidth ||
            neighborCellY < 0 || neighborCellY >= gridHeight)
        {
          continue;
        }

        int cellIdx = getCellIndex(neighborCellX, neighborCellY);
        auto it = cells.find(cellIdx);
        if (it != cells.end())
        {
          // Filter particles by actual distance
          for (sph::Particle *p : it->second)
          {
            float dx = x - p->getPosition().x;
            float dy = y - p->getPosition().y;
            float distSquared = dx * dx + dy * dy;

            if (distSquared <= radiusSquared)
            {
              neighbors.push_back(p);
            }
          }
        }
      }
    }

    return neighbors;
  }

  std::vector<sph::Particle *> Grid::getNeighbors(const sph::Particle *particle,
                                                  float radius)
  {
    auto pos = particle->getPosition();
    return getNeighbors(pos.x, pos.y, radius);
  }

  int Grid::getCellIndex(float x, float y) const
  {
    auto [cellX, cellY] = positionToCell(x, y);
    return getCellIndex(cellX, cellY);
  }

  int Grid::getCellIndex(int cellX, int cellY) const
  {
    return cellY * gridWidth + cellX;
  }

  std::pair<int, int> Grid::positionToCell(float x, float y) const
  {
    int cellX = static_cast<int>(x * invCellSize);
    int cellY = static_cast<int>(y * invCellSize);
    return {cellX, cellY};
  }

} // namespace sph
