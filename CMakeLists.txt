cmake_minimum_required(VERSION 3.15)
project(DensityMatrixQuantumSimulator)

option(DMQS_TESTS "DMQSUnit Test" TRUE)

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_CXX_STANDARD_EXTENSIONS FALSE)

find_package(Armadillo REQUIRED)
find_package(spdlog REQUIRED)

add_library(external SHARED gates.cpp dmqs.cpp)
target_compile_features(external PRIVATE cxx_std_23)
target_compile_options(external PUBLIC -Wpedantic -Wall -Wextra)
target_link_libraries(external PUBLIC armadillo spdlog)


if (DMQS_TESTS)
  add_executable(dmqs_test dmqs_test.cpp)
  target_link_libraries(dmqs_test external)
  add_custom_command(
    TARGET dmqs_test
    POST_BUILD
    COMMAND dmqs_test
  )
endif(DMQS_TESTS)