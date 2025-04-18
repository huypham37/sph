cmake_minimum_required(VERSION 3.10)
project(SPH_Simulation VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configure build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Find SFML 3.0
find_package(SFML 3 COMPONENTS Graphics Window System REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add source files
file(GLOB SOURCES "src/*.cpp")

# Create executable
add_executable(sph_simulation ${SOURCES})

# Link libraries using new SFML 3.0 target names
target_link_libraries(sph_simulation SFML::Graphics SFML::Window SFML::System)

# Install target
install(TARGETS sph_simulation DESTINATION bin)