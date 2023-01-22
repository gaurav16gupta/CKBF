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
    uint32_t N = 10;
    string filellistname = config.fastqFileName;
    uint64_t MAXRANGE = 8589934592; //2^33
    bool insert =0;
    bool query =1;
    // Create and open a log file
    // ofstream InsertLogs("./results/insertLogs.txt");
    Hasher *hasher[config.numThreads]; // each thread gets its own hasher
        for (uint32_t i = 0; i < config.numThreads; ++i) {
            hasher[i] = config.hashType == Config::MURMUR_HASH
            ? static_cast<Hasher*>(new MurmurHasher(MAXRANGE, config.k, config.seed))
            : (config.hashType == Config::BRUTE_FORCE_FUZZY_HASH
                ? static_cast<Hasher*>(new FuzzyHasher(MAXRANGE, config.k, config.kMer, config.universalHashRange, config.seed))
                : static_cast<Hasher*>(new EfficientFuzzyHasher(MAXRANGE, config.k, config.kMer, config.universalHashRange, config.seed)));
        }

    BloomFilter** ArBF_array; 
    ArBF_array = new BloomFilter*[N]; //array of pointers
    vector<uint32_t> numInsert;
    vector<uint32_t> ranges;
    uint32_t range;
    uint32_t n = 0;

    ifstream filellist(filellistname);
    string fileName;

    if (insert){
        while (getline(filellist, fileName) && n<=N ) {
            vector<string> sequences = getFastqData("./data/fastqFiles/"+ fileName);
            range = config.rangefactor*config.k*sequences.size();
            ArBF_array[n] = new BloomFilter(range, config.k, config.disk, fileName.substr(0,fileName.length() - 6)+ " W");
            printf("%s, #Inserts:%ld, BF size:%d", fileName.c_str(), sequences.size(), range);
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
                    ArBF_array[n]->insert(hashes);
                    // t3 = chrono::high_resolution_clock::now();
                    // hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
                    // bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
                    // counter[n] += 1;
                }
            }
            ArBF_array[n]->release();
            n++;
            numInsert.push_back(sequences.size());
            ranges.push_back(range);
        }
        if (n<N) printf("WARNING: less than %d .fastq files!", N);
        write_vector("./results/ranges.txt", ranges);
    }

    if(query){
        string queryFile = config.queryFileName;
        uint32_t nq = 1000;
        vector<string> queries(nq);
        vector<vector<uint32_t>> GT(nq);
        getQueryforArBF(queryFile, queries, GT, N);

        vector<vector<uint32_t>> q_results(nq);
        uint32_t fpCount,fpc = 0;
        uint32_t hashTimeAccu = 0, bfTimeAccu = 0;
        chrono::time_point<chrono::high_resolution_clock> t1, t2, t3;

        ranges= read_vector("./results/ranges.txt");
        while (getline(filellist, fileName) && n<=N ) {
            range = ranges[n];
            ArBF_array[n] = new BloomFilter(range, config.k, config.disk, fileName.substr(0,fileName.length() - 6)+ " R");
            n++;
            cout<<fileName.substr(0,fileName.length() - 6)+ " R"<<endl;
        }
        // # pragma omp parallel for 
        for (size_t i = 0; i < nq; ++i) {
            uint32_t threadId = 0; // omp_get_thread_num();
            uint64_t hashes[config.k];
            hasher[threadId]->setSequence(queries[i]);
            while (hasher[threadId]->hasNext()) {
                t1 = chrono::high_resolution_clock::now();
                hasher[threadId]->hash(hashes);
                t2 = chrono::high_resolution_clock::now();
                for (size_t n = 0; n < N; ++n){
                    for (uint8_t k = 0; k < config.k; ++k){
                            hashes[k] = hashes[k]%ranges[n];
                    }
                    if (ArBF_array[n]->test(hashes)){
                        q_results[i].push_back(n);
                    }
                }
                cout<<"3"<<endl;
                t3 = chrono::high_resolution_clock::now();
                hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
                bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
            }
            cout<<"4"<<endl;
            fpc = (q_results[i].size()- GT[i].size());
            assert (fpc>=0);
            fpCount += fpc;
        }
        fpCount = fpCount/nq;

        cout << "Hash Time: " << hashTimeAccu << "; Read Time: " << bfTimeAccu << endl;
        cout << "Average False Positive Rate: " << static_cast<float>(fpCount) / N << endl;
    }
    return 0;
}

//todo: 
// 1) choice of ABF experiment in RAM, use disk flag (if 0) to save and load BF arrays n RAM 