cmake_minimum_required(VERSION 3.22.1)
project(TSP C)

set(CMAKE_C_STANDARD 99)

# Find OpenMP package
find_package(OpenMP REQUIRED)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build_releas)

set(SOURCE_FILES tsp.c queue.c queue.h tscp.h)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
        # Add the two source files
add_executable(TSP ${SOURCE_FILES})

# Link with OpenMP
target_link_libraries(TSP PRIVATE OpenMP::OpenMP_C)

#O3 flag
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")