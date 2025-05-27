#!/bin/bash
# Build the SPH project
source ./setup_env.sh

cd "/home/coder/project/sph/build"
cmake -DCMAKE_BUILD_TYPE=Release       -DSFML_DIR="${SFML_DIR}/lib/cmake/SFML"       -DCMAKE_CXX_FLAGS="-O3 -march=native -fopenmp"       ..

make -j${OMP_NUM_THREADS}

echo "Build complete. You can run the application with:"
echo "cd /home/coder/project/sph/build && ./sph_simulation"
