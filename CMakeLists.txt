cmake_minimum_required(VERSION 3.6)
project(rds)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
file(GLOB_RECURSE HEADER_FILES src/*.h)
add_executable(rds ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(rds PUBLIC src)
