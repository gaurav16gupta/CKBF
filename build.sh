#!/bin/bash

g++ src/MurmurHash3.cpp src/BloomFilter.cpp src/main.cpp -Wall -fopenmp -std=c++17 -O3 -Iinclude/ -g -o debug/program
g++ src/MurmurHash3.cpp src/BloomFilter.cpp src/ArBF.cpp -Wall -fopenmp -std=c++17 -O3 -Iinclude/ -g -o debug/abf