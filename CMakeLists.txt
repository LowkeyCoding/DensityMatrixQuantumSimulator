cmake_minimum_required(VERSION 3.15)
project(DensityMatrixQuantumSimulator)

option(GATES_TESTS "GATESUnit Test" FALSE)
option(DMQS_TESTS "DMQSUnit Test" TRUE)
option(UPPAAL_TESTS "UPPAALUnit Test" TRUE)

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_CXX_STANDARD_EXTENSIONS FALSE)

# Global compile options
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic -Wall -Wextra")

find_package(Armadillo REQUIRED)
find_package(spdlog REQUIRED)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(bindings)
add_subdirectory(tests)