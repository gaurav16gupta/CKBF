#!/bin/bash

g++ src/MurmurHash3.cpp src/BloomFilter.cpp src/main.cpp -Wall -fopenmp -std=c++17 -O2 -Iinclude/ -g -o debug/program