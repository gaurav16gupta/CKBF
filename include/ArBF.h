#ifndef _ArBF_
#define _ArBF_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
// #include <omp.h>
// #include "BloomFilter.h"
// #include "utils.h"
// #include "config.h"
// #include "Hasher.h"
#include <assert.h>
#include <cstdlib>
#include <dirent.h>

using namespace std;

class ArBF{
public:
    ArBF(uint32_t N);
    void insert(string folder, const Config config);
    BloomFilter** ArBF_array; 
    uint32_t N;
    uint32_t *numInsert;
    uint32_t *rangeOffset;
};

#endif

