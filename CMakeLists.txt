cmake_minimum_required(VERSION 3.24)
project(TSP C)

set(CMAKE_C_STANDARD 99)


# Find OpenMP package
find_package(OpenMP REQUIRED)

set(SOURCE_FILES tsp.c queue.c queue.h tscp.h)
# Add the two source files
add_executable(TSP tsp.c queue.c queue.h tscp.h)

# Link with OpenMP
target_link_libraries(TSP PRIVATE OpenMP::OpenMP_C)
# Fix undefined reference of functoins in the queue.c file


# Fix undefined reference to `omp_get_wtime`
#target_link_libraries(TSP PRIVATE m)

# Set compiler flags
target_compile_options(TSP PRIVATE -O3)