cmake_minimum_required(VERSION 3.5)
project(rds)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp mesh.cpp point2D.cpp vector2D.cpp)
add_executable(rds ${SOURCE_FILES})