cmake_minimum_required(VERSION 3.24.0)
project(tsp C )

set(CMAKE_C_STANDARD 99)

# Find OpenMP package
find_package(OpenMP REQUIRED)

#set(CMAKE_EXE_LINKER_FLAGS "-static")

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build_releas)

set(SOURCE_FILES tsp.c queue.c queue.h tscp.h)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
        # Add the two source files
add_executable(tsp ${SOURCE_FILES})

# Link with OpenMP
target_link_libraries(tsp PRIVATE OpenMP::OpenMP_C)

#O3 flag
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -g  ")