cmake_minimum_required(VERSION 3.16)
project(idl)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -fopenmp -O2 -g")

set(SOURCE_FILES src/main.cpp src/MurmurHash3.cpp src/BloomFilter.cpp)

add_executable(single_bf ${SOURCE_FILES})
target_include_directories(single_bf PUBLIC include)
