#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
// #include <omp.h>
#include "BloomFilter.h"
#include "utils.h"
#include "config.h"
#include "Hasher.h"
#include <filesystem>

using namespace std;

int main(int argc, char* argv[]) {
    const Config config = getConfigs("./configs/" + (argc >= 2 ? string(argv[1]) : "default") + ".cfg");
    config.print();

    // std::string fastqFileName = config.fastqFileName;
    // std::string queryFileName = config.queryFileName;
    // if(argc >= 2){fastqFileName = argv[1];}
    // if(argc >= 3){queryFileName = fastqFileName+ "query.p"+ argv[2];}
    // fastqFileName = "data/fastqFiles/" + fastqFileName + ".fastq";
    // queryFileName = "data/queries/" + queryFileName;
    // vector<string> sequences = getFastqData(fastqFileName);
    // vector<string> querySequences = getQueryData(queryFileName);
    // uint32_t range = 10*(203-31)*filesystem::file_size(fastqFileName)/(203*2 + 26*2);
    // // range = (range/(4096*8));
    // // range = range*(4096*8);
    // cout << "range: "<<range<<" bits, or "<<range/(8*1024*1024)<<" Mbs"<<endl;
    // BloomFilter bf(range, config.k, 1);

    vector<string> sequences = getFastqData("../data/" + config.fastqFileName + ".fastq");
    // vector<string> querySequences = getQueryData("../data/" + config.queryFileName);
    BloomFilter bf(config.range, config.k, 1);

    // omp_set_num_threads(config.numThreads);

    Hasher *hasher[config.numThreads]; // each thread gets its own hasher
    for (uint32_t i = 0; i < config.numThreads; ++i) {
        hasher[i] = config.hashType == Config::MURMUR_HASH
        ? static_cast<Hasher*>(new MurmurHasher(static_cast<uint32_t>(config.range), config.k, config.seed))
        : static_cast<Hasher*>(new EfficientFuzzyHasher(config.range, config.k, config.kMer, config.universalHashRange, config.seed));
    }
    // uint32_t hashTimeAccu = 0, bfTimeAccu = 0;
    chrono::time_point<chrono::high_resolution_clock> t1, t2, t3;
    t1 = chrono::high_resolution_clock::now();
    // # pragma omp parallel for 
    for (size_t i = 0; i < sequences.size(); ++i) {
        uint32_t hashes[config.k];
        uint32_t threadId = 0; // omp_get_thread_num();
        hasher[threadId]->setSequence(sequences[i]);
        while (hasher[threadId]->hasNext()) {
            hasher[threadId]->hash(hashes);
            bf.insert(hashes);
            // t3 = chrono::high_resolution_clock::now();
            // hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
            // bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
        }
    }
    t2 = chrono::high_resolution_clock::now();
    cout << "Insert Time: " << chrono::duration_cast<chrono::microseconds>(t2 - t1).count()/1000000.0 << endl;

    // hashTimeAccu = 0, bfTimeAccu = 0;
    uint32_t fpCount = 0;
    t1 = chrono::high_resolution_clock::now();
    // # pragma omp parallel for 
    for (size_t i = 0; i < sequences.size(); ++i) {
        uint32_t threadId = 0; // omp_get_thread_num();
        uint32_t hashes[config.k];
        sequences[i][10] = 'B';
        hasher[threadId]->setSequence(sequences[i]);
        chrono::time_point<chrono::high_resolution_clock> t1, t2, t3;
        bool fp = true;
        while (hasher[threadId]->hasNext()) {
            hasher[threadId]->hash(hashes);
            // t2 = chrono::high_resolution_clock::now();
            fp = fp && bf.test(hashes);
            // t3 = chrono::high_resolution_clock::now();
            // hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
            // bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
        }
        fpCount += (fp ? 1 : 0);
    }
    t2 = chrono::high_resolution_clock::now();
    cout << "Query Time: " << chrono::duration_cast<chrono::microseconds>(t2 - t1).count()/1000000.0 << endl;
    cout << "False Positive Rate: " << static_cast<float>(fpCount) / sequences.size() << endl;
    // bf.release();
    return 0;
}
