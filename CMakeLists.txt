cmake_minimum_required(VERSION 3.15)
project(DensityMatrixQuantumSimulator)

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_CXX_STANDARD_EXTENSIONS FALSE)
set(DMQS_EX TRUE)

# Global compile options
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic -Wall -Wextra -O3")  # use -L/opt/cuda/lib64 -lnvblas for gpu support 

include(cmake/sanitizers.cmake)
include(cmake/armadillo.cmake)
find_package(OpenBLAS REQUIRED)

include_directories(include)
# Add subdirectories
add_subdirectory(src)
add_subdirectory(bindings)

include(cmake/doctest.cmake)
enable_testing()
add_subdirectory(tests)


if(DMQS_EX)
    add_executable(densecoding examples/densecoding.cpp)
    target_link_libraries(densecoding dmqs_uppaal)
endif()