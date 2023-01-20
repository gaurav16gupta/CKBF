#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include "BloomFilter.h"
#include "utils.h"
#include "config.h"
#include "Hasher.h"
#include <assert.h>
#include <cstdlib>
#include <dirent.h>

using namespace std;

int main(int argc, char** argv) {
    const Config config = getConfigs(argv[1]);
    config.print();
    uint32_t N;
    string filellistname;

    uint32_t numInsert[N];
    uint32_t rangeOffset[N];
    BloomFilter** ArBF_array; 
    ArBF_array = new BloomFilter*[N]; //array of pointers

    uint32_t range;
    uint32_t n = 0;

    ifstream filellist(filellistname);
    string fileName;
    while (getline(filellist, fileName) && n<=N ) {
        vector<string> sequences = getFastqData("./data/"+ fileName);
        range = config.rangefactor*config.k*sequences.size();
        ArBF_array[n] = new BloomFilter(range, config.k, config.disk, fileName.substr(0,fileName.length() - 6));

        Hasher *hasher[config.numThreads]; // each thread gets its own hasher
        for (uint32_t i = 0; i < config.numThreads; ++i) {
            hasher[i] = config.hashType == Config::MURMUR_HASH
            ? static_cast<Hasher*>(new MurmurHasher(static_cast<uint32_t>(config.range), config.k, config.seed))
            : (config.hashType == Config::BRUTE_FORCE_FUZZY_HASH
                ? static_cast<Hasher*>(new FuzzyHasher(static_cast<uint32_t>(config.range), config.k, config.kMer, config.universalHashRange, config.seed))
                : static_cast<Hasher*>(new EfficientFuzzyHasher(static_cast<uint32_t>(config.range), config.k, config.kMer, config.universalHashRange, config.seed)));
        }
        // # pragma omp parallel for 
        for (size_t i = 0; i < sequences.size(); ++i) {
            uint32_t hashes[config.k];
            uint32_t threadId = 0; // omp_get_thread_num();
            hasher[threadId]->setSequence(sequences[i]);
            while (hasher[threadId]->hasNext()) {
                // t1 = chrono::high_resolution_clock::now();
                hasher[threadId]->hash(hashes);
                // t2 = chrono::high_resolution_clock::now();
                ArBF_array[n]->insert(hashes);
                // t3 = chrono::high_resolution_clock::now();
                // hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
                // bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
                // counter[n] += 1;
            }
        }
        n++;
        numInsert[n] = sequences.size();
        rangeOffset[n+1] = rangeOffset[n] + range;
        printf("%s\t", fileName.c_str());
    }
    if (n<N) printf("WARNING: less than %d .fastq files!", N);

    return 0;
}