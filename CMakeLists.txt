cmake_minimum_required(VERSION 3.6)
project(rds)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -ipo -xSSE4.2")

file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
file(GLOB_RECURSE HEADER_FILES src/*.h)
add_executable(rds ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(rds PUBLIC src)

set_property(TARGET rds PROPERTY CXX_STANDARD 14)
set_property(TARGET rds PROPERTY CXX_STANDARD_REQUIRED ON)
