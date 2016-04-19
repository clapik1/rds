cmake_minimum_required(VERSION 3.5)
project(rds)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp mesh.cpp mesh.h geo/point2D.cpp geo/point2D.h geo/vector2D.cpp geo/vector2D.h geo/triangle2D.cpp geo/triangle2D.h solver.cpp solver.h)
add_executable(rds ${SOURCE_FILES})