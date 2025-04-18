- adding threadpool. 
- parallel like using thread pool
	- How current implementation works with 


naturally parallizable 
Density and pressure calculation (computeDensityPressure())
Force computation (computeForces())
Position/velocity integration (integrate())
Collision resolution (resolveCollisions())


Using std::unordered_map vs. Vector of Vectors:

We're currently using an unordered_map for cell storage, which has lookup overhead.
A 2D array or vector of vectors might be more efficient if the grid dimensions are reasonable.
SIMD Optimizations:

For distance calculations, SIMD instructions could compute multiple particle distances in parallel.
More Efficient Collision Detection:

The current approach checks particle collisions after the grid has been updated for rendering.
We could integrate collision detection directly within the SPH force calculation loop.
