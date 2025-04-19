- Verlet list (sec 2)


Paralellizable area: 
	naturally parallizable 
	- Density and pressure calculation (computeDensityPressure())
	Force computation (computeForces())
	Position/velocity integration (integrate())
	Collision resolution (resolveCollisions())

Bottlenecks that can improved

	Using std::unordered_map vs. Vector of Vectors:

	We're currently using an unordered_map for cell storage, which has lookup overhead.
	A 2D array or vector of vectors might be more efficient if the grid dimensions are reasonable.
	SIMD Optimizations:

	For distance calculations, SIMD instructions could compute multiple particle distances in parallel.
	More Efficient Collision Detection:

	The current approach checks particle collisions after the grid has been updated for rendering.
	We could integrate collision detection directly within the SPH force calculation loop.

#file:.system_prompt.md  the getNeighbour is the expensive task, why is it called in many places? 
i think it should only be called one in each step of the simulation right? and how about paralllel that tasks? we already have domain decomposer which decompose the program into smaller subdomain and subdoamin has particles, why dont we using different threads to acceelarte searching