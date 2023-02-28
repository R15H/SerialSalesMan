cmake_minimum_required(VERSION 3.24)
project(TSP C)

set(CMAKE_C_STANDARD 99)

add_executable(untitled tsp.c tscp.h)

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Add the two source files
add_executable(TSP tsp.c queue.c queue.h tscp.h)

# Link with OpenMP
target_link_libraries(TSP PRIVATE OpenMP::OpenMP_C)

# Set compiler flags
target_compile_options(TSP PRIVATE -O3)