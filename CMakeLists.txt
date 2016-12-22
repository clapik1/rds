cmake_minimum_required(VERSION 3.6)
project(rds)
#find_package(CUDA)

#list(APPEND CUDA_NVCC_FLAGS "-arch=sm_20;-std=c++11;-O2;-DVERBOSE")
#list(APPEND CMAKE_CXX_FLAGS "-std=c++14 -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

#set(CUDA_PROPAGATE_HOST_FLAGS OFF)

#set(CUDA_COMMON_DIR /home/michal/workspace/cuda_samples/NVIDIA_CUDA-8.0_Samples/common/inc)
#include_directories(${CUDA_COMMON_DIR})
#set(HEADER_FILES ${CUDA_COMMON_DIR}/helper_cuda.h)

file(GLOB_RECURSE SOURCE_FILES src/*.cpp) #src/*.cu
#set(SOURCE_FILES src/*.cpp)
add_executable(rds ${SOURCE_FILES})
