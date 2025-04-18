cmake_minimum_required(VERSION 3.10)
project(SPH_Simulation VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configure build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# OpenMP configuration for macOS/Apple Clang
if(APPLE AND CMAKE_C_COMPILER_ID MATCHES "Clang")
    # Check if OpenMP is installed via Homebrew
    if(EXISTS "/opt/homebrew/opt/libomp")
        # OpenMP installed via Homebrew on Apple Silicon Mac
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    elseif(EXISTS "/usr/local/opt/libomp")
        # OpenMP installed via Homebrew on Intel Mac
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
    else()
        message(WARNING "OpenMP not found via Homebrew. Please install with: brew install libomp")
    endif()
endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Find SFML 3.0
find_package(SFML 3 COMPONENTS Graphics Window System REQUIRED)

# Find OpenMP with additional variables now set
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - Version: ${OpenMP_CXX_VERSION}")
    message(STATUS "OpenMP C++ flags: ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP C++ libraries: ${OpenMP_CXX_LIBRARIES}")
else()
    message(FATAL_ERROR "OpenMP not found but required - aborting build")
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add source files
file(GLOB SOURCES "src/*.cpp")

# Create executable
add_executable(sph_simulation ${SOURCES})

# Link libraries using new SFML 3.0 target names
target_link_libraries(sph_simulation SFML::Graphics SFML::Window SFML::System OpenMP::OpenMP_CXX)

# Install target
install(TARGETS sph_simulation DESTINATION bin)