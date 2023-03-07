cmake_minimum_required(VERSION 3.22)
project(TSP C)

set(CMAKE_C_STANDARD 99)

# Find OpenMP package
find_package(OpenMP REQUIRED)

set(SOURCE_FILES tsp.c queue.c queue.h tscp.h)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

        # Add the two source files
add_executable(TSP ${SOURCE_FILES})

# Link with OpenMP
target_link_libraries(TSP PRIVATE OpenMP::OpenMP_C)

target_comPile_options(TSP PRIVATE -O3)