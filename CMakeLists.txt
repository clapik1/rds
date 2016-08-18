cmake_minimum_required(VERSION 3.5)
project(rds)

# enable fortran, c, and c++ language
enable_language(CXX)

# project version
set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_PATCH 0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
#set(SOURCE_FILES src/*.cpp)
add_executable(rds ${SOURCE_FILES})
