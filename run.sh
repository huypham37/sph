#!/bin/bash
# filepath: run_sph_cluster.sh

# Exit on any error
set -e

# ===== Configuration =====
# Adjust these according to your cluster environment
NUM_CORES=24            # Number of cores to use for OpenMP
CLUSTER_QUEUE="compute"  # Your cluster's job queue
JOB_NAME="sph_simulation" # Name for the job
WALLTIME="01:00:00"     # Maximum run time (HH:MM:SS)

# ===== Directory Setup =====
WORK_DIR=$(pwd)
BUILD_DIR="${WORK_DIR}/build"
SOURCE_DIR="${WORK_DIR}"
DEPS_DIR="${WORK_DIR}/deps"

mkdir -p "${BUILD_DIR}"
mkdir -p "${DEPS_DIR}"

echo "===== Starting SPH Simulation Cluster Job ====="
echo "Working directory: ${WORK_DIR}"
echo "Building in: ${BUILD_DIR}"

# ===== Dependencies =====
# SFML is required for the graphics
echo "Installing dependencies..."

# Install SFML from source if needed
if [ ! -d "${DEPS_DIR}/SFML-install" ]; then
    cd "${DEPS_DIR}"
    echo "Downloading SFML..."
    wget https://github.com/SFML/SFML/archive/refs/tags/2.5.1.tar.gz -O SFML-2.5.1.tar.gz
    tar -xzf SFML-2.5.1.tar.gz
    cd SFML-2.5.1
    
    # Build SFML
    echo "Building SFML..."
    mkdir -p build && cd build
    cmake -DCMAKE_INSTALL_PREFIX="${DEPS_DIR}/SFML-install" -DBUILD_SHARED_LIBS=ON ..
    make -j${NUM_CORES}
    make install
    
    echo "SFML installed successfully."
    cd "${WORK_DIR}"
fi

# ===== Configure OpenMP =====
# Set OpenMP environment variables
export OMP_NUM_THREADS=${NUM_CORES}
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_SCHEDULE="dynamic,64"

# ===== Build Project =====
echo "Building SPH simulation..."
cd "${BUILD_DIR}"

# Generate build files with CMake
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH="${DEPS_DIR}/SFML-install" \
      -DCMAKE_CXX_FLAGS="-O3 -march=native -fopenmp" \
      "${SOURCE_DIR}"

# Compile
make -j${NUM_CORES}

echo "Build completed successfully."

# ===== Run Simulation =====
echo "Starting SPH simulation with ${NUM_CORES} threads..."

# Optional: Create job submission script for cluster systems like SLURM
cat > run_job.sh << EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${NUM_CORES}
#SBATCH --time=${WALLTIME}
#SBATCH --partition=${CLUSTER_QUEUE}

cd ${BUILD_DIR}

# Set OpenMP environment
export OMP_NUM_THREADS=${NUM_CORES}
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_SCHEDULE="dynamic,64"

# Run with performance measurements
time ./sph_simulation

echo "Simulation completed."
EOF

chmod +x run_job.sh

echo "===== Job submission script created ====="
echo "To submit the job:"
echo "  sbatch run_job.sh"
echo ""
echo "To run directly:"
echo "  ./sph_simulation"
echo ""
echo "Setup complete!"