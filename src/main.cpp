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
#include <assert.h>
#include <cstdlib>

using namespace std;


string poison(int str_id, string s, int strength) {
    srand(10256*str_id);
    for(int i=0; i < strength; i++) {
        s[rand()%s.length()] = 'X';
    }
    return s;
}

int main(int argc, char** argv) {
    const Config config = getConfigs(argv[1]);
    config.print();

    vector<string> sequences = getFastqData("/home/apd10/code/CKBF/data/" + config.fastqFileName + ".fastq");
    // vector<string> querySequences = getQueryData("../data/" + config.queryFileName);
    BloomFilter bf(config.range, config.k);

    // omp_set_num_threads(config.numThreads);
    assert(config.numThreads == 1);
    Hasher *hasher[config.numThreads]; // each thread gets its own hasher
    for (uint32_t i = 0; i < config.numThreads; ++i) {
        hasher[i] = config.hashType == Config::MURMUR_HASH
        ? static_cast<Hasher*>(new MurmurHasher(static_cast<uint32_t>(config.range), config.k, config.seed))
        : static_cast<Hasher*>(new EfficientFuzzyHasher(static_cast<uint32_t>(config.range), config.k, config.kMer, config.universalHashRange, config.seed));
    }

    uint32_t hashTimeAccu = 0, bfTimeAccu = 0, counter = 0;
    // # pragma omp parallel for reduction(+:hashTimeAccu, bfTimeAccu, counter)
    for (size_t i = 0; i < sequences.size(); ++i) {
        uint32_t hashes[config.k];
        uint32_t threadId = 0; // omp_get_thread_num();
        hasher[threadId]->setSequence(sequences[i]);
        chrono::time_point<chrono::high_resolution_clock> t1, t2, t3;
        while (hasher[threadId]->hasNext()) {
            t1 = chrono::high_resolution_clock::now();
            hasher[threadId]->hash(hashes);
            t2 = chrono::high_resolution_clock::now();
            bf.insert(hashes);
            t3 = chrono::high_resolution_clock::now();
            hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
            bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
            counter += 1;
        }
    }

    cout << "Hash Time: " << hashTimeAccu << "; BF Insert Time: " << bfTimeAccu << endl;
    cout << "BF Packing: " << bf.count() << '/' << bf.size << endl;
    cout << "# of items inserted: " << counter << endl;

    hashTimeAccu = 0, bfTimeAccu = 0;
    uint32_t fpCount = 0;
    // # pragma omp parallel for reduction(+:hashTimeAccu, bfTimeAccu, fpCount)
    for (size_t i = 0; i < sequences.size(); ++i) {
        uint32_t threadId = 0; // omp_get_thread_num();
        uint32_t hashes[config.k];
        sequences[i] = poison(i, sequences[i], config.poison);
        hasher[threadId]->setSequence(sequences[i]);
        chrono::time_point<chrono::high_resolution_clock> t1, t2, t3;
        bool fp = true;
        while (hasher[threadId]->hasNext()) {
            t1 = chrono::high_resolution_clock::now();
            hasher[threadId]->hash(hashes);
            t2 = chrono::high_resolution_clock::now();
            fp = fp && bf.test(hashes);
            t3 = chrono::high_resolution_clock::now();
            hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
            bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
        }
        fpCount += (fp ? 1 : 0);
    }
    cout << "Hash Time: " << hashTimeAccu << "; BF Query Time: " << bfTimeAccu << endl;
    cout << "Positive Rate: " << static_cast<float>(fpCount) / sequences.size() << endl;
    return 0;
}
