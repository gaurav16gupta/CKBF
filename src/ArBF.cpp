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
#include <stdio.h>

using namespace std;

int main(int argc, char** argv) {
    const Config config = getConfigs(argv[1]);
    config.print();
    string filellistname = config.fastqFileName;
    bool insert = !config.queryOnly;
    bool query = true;
    // Create and open a log file
    // ofstream InsertLogs("./results/insertLogs.txt");
    Hasher *hasher[config.numThreads]; // each thread gets its own hasher
    for (uint32_t i = 0; i < config.numThreads; ++i) {
        hasher[i] = config.hashType == Config::MURMUR_HASH
        ? static_cast<Hasher*>(new MurmurHasher(UINT64_MAX, config.k, config.seed))
        : (config.hashType == Config::BRUTE_FORCE_FUZZY_HASH
            ? static_cast<Hasher*>(new FuzzyHasher(UINT64_MAX, config.k, config.kMer, config.universalHashRange, config.seed))
            : static_cast<Hasher*>(new EfficientFuzzyHasher(UINT64_MAX, config.k, config.kMer, config.universalHashRange, config.seed)));
    }

    BloomFilter** ArBF_array; 
    ArBF_array = new BloomFilter*[config.firstNFiles]; //array of pointers
    vector<uint32_t> numInsert;
    vector<uint32_t> ranges;
    uint32_t range;
    uint32_t numFiles = 0;

    string fileName;

    if (insert) {
        ifstream filellist(filellistname);
        while (numFiles < config.firstNFiles && getline(filellist, fileName)) {
            vector<string> sequences = getFastqData("/scratch1/gg29/CKBF/data/fastqFiles/" + fileName);
            range = config.rangefactor*config.k*sequences.size();
            ArBF_array[numFiles] = new BloomFilter(range, config.k, config.disk, fileName.substr(0,fileName.length() - 6)+ " W");
            // # pragma omp parallel for 
            for (size_t i = 0; i < sequences.size(); ++i) {
                uint64_t hashes[config.k];
                uint32_t threadId = 0; // omp_get_thread_num();
                hasher[threadId]->setSequence(sequences[i]);
                while (hasher[threadId]->hasNext()) {
                    // t1 = chrono::high_resolution_clock::now();
                    hasher[threadId]->hash(hashes);
                    for (uint8_t k = 0; k < config.k; ++k){
                        hashes[k] = hashes[k]%range;
                    }
                    
                    // t2 = chrono::high_resolution_clock::now();
                    ArBF_array[numFiles]->insert(hashes);
                    // t3 = chrono::high_resolution_clock::now();
                    // hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
                    // bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
                    // counter[n] += 1;
                }
            }
            cout << "BF #" << numFiles << " packing: " << ArBF_array[numFiles]->count() << '/' << ArBF_array[numFiles]->size << "; # of inserts: " << sequences.size() << endl;
            ArBF_array[numFiles]->release();
            numFiles++;
            numInsert.push_back(sequences.size());
            ranges.push_back(range);
        }
        cout << "Indexed " << numFiles << " .fastq files." << endl;
        write_vector("./results/ranges.txt", ranges);
    }

    uint32_t n = 0;
    if (query) {
        string queryFile = config.queryFileName;
        uint32_t nq = 1000; // TODO: remove hardcoded nq
        vector<string> queries(nq);
        vector<vector<uint32_t>> GT(nq);
        getQueryforArBF(queryFile, queries, GT, numFiles);

        uint32_t hashTimeAccu = 0, bfTimeAccu = 0;
        chrono::time_point<chrono::high_resolution_clock> t1, t2;

        ranges= read_vector("./results/ranges.txt");
        ifstream filellist(filellistname);
        while (n < numFiles && getline(filellist, fileName)) {
            range = ranges[n];
            ArBF_array[n] = new BloomFilter(range, config.k, config.disk, fileName.substr(0,fileName.length() - 6)+ " R");
            n++;
        }

        ofstream queryResults("./results/queryResults.txt", ofstream::trunc);

        // # pragma omp parallel for 
        for (size_t i = 0; i < nq; ++i) {
            uint32_t threadId = 0; // omp_get_thread_num();
            hasher[threadId]->setSequence(queries[i]);
            uint32_t nHashes = hasher[threadId]->nHashesInSequence();
            uint64_t hashes[config.k * nHashes];
            
            t1 = chrono::high_resolution_clock::now();
            for (uint32_t l = 0; l < nHashes; ++l) {
                hasher[threadId]->hash(hashes + l * config.k);
            }
            t2 = chrono::high_resolution_clock::now();
            hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();

            queryResults << queries[i];

            for (size_t n = 0; n < numFiles; ++n){
                bool positive = true;
                t1 = chrono::high_resolution_clock::now();
                for (uint32_t l = 0; l < nHashes; ++l) {
                    positive = positive && ArBF_array[n]->test(hashes, ranges[n]);
                }
                t2 = chrono::high_resolution_clock::now();
                bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
                if (positive) {
                    queryResults << ',' << n;
                }
            }
            queryResults << endl;
        }
        cout << "Hash Time: " << hashTimeAccu << "; BF Lookup Time: " << bfTimeAccu << endl;

        for (uint32_t i = 0; i < numFiles; ++i) {
            ArBF_array[i]->release();
        }
    }
    return 0;
}

//todo: 
// 1) choice of ABF experiment in RAM, use disk flag (if 0) to save and load BF arrays n RAM 
// 2) stop when first kmer gives False, or calculate the %  kmer content
// 3) get query,fileID queryfile