cmake_minimum_required(VERSION 3.12)
project(libpynamd LANGUAGES CXX)

# Compiler options: C++17, warnings, and optimization
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O3)              # Generic maximum code optimization
add_compile_options(-ftree-vectorize) # SIMD vectorize the loops
add_compile_options(-funroll-loops)   # explicit optimize the loops
add_compile_options(-Wall -Wextra -Wpedantic) # Enable warning
# add archetecture specific optimization
add_compile_options(-march=native)

# add_include direcotry
include_directories(${CMAKE_SOURCE_DIR}/include)


# Find Eigen3 and pybind11
set(EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/eigen/)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/pybind11)


# Add binding library
add_subdirectory(libpynamd)
