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

char rep(char x) {
    return 'X';
    // Make a sure shot false positive
    switch(x) {
        case 'A':
            return 'C';
        case 'C':
            return 'G';
        case 'G':
            return 'T';
        case 'T':
            return 'A';
        default:
            return 'X';
    }
}

string poison(int str_id, string s, int strength) {
    srand(10256*str_id);
    int l = max((size_t)32, rand() % s.length());
    int start = min((size_t)(l-32), rand() % s.length());
    s = s.substr(start, l);
    int temp;
    for(int i=0; i < strength; i++) {
        temp = rand()%s.length();
        s[temp] = rep(s[temp]);
    }
    return s;
}

int main(int argc, char** argv) {
    const Config config = getConfigs(argv[1], argc - 2, argv + 2);
    config.print();

    vector<string> sequences = getFastqData("./data/fastqFiles/" + config.fastqFileName + ".fastq");
    // vector<string> querySequences = getQueryData("../data/" + config.queryFileName);
    BloomFilter bf(config.range, config.k, config.disk, config.fastqFileName+ " W");
    // omp_set_num_threads(config.numThreads);
    assert(config.numThreads == 1);
    Hasher *hasher[config.numThreads]; // each thread gets its own hasher
    for (uint32_t i = 0; i < config.numThreads; ++i) {
        switch(config.hashType) {
            case Config::MURMUR_HASH:
                hasher[i] = static_cast<Hasher*>(new MurmurHasher(config.range, config.k, config.seed));
                break;
            case Config::BRUTE_FORCE_FUZZY_HASH:
                hasher[i] = static_cast<Hasher*>(new FuzzyHasher(config.range, config.k, config.kMer, config.universalHashRange, config.seed));
                break;
            case Config::FUZZY_HASH:
                hasher[i] = static_cast<Hasher*>(new EfficientFuzzyHasher(config.range, config.k, config.kMer, config.universalHashRange, config.seed));
                break;
            case Config::FUZZY_HASH_EXP:
                hasher[i] = static_cast<Hasher*>(new EfficientFuzzyHasherEXP(config.range, config.k, config.kMer, config.universalHashRange, config.seed));
                break;
            default:
                cout << "HASH TYPE not recognized" << endl;
                return 0;
        }
    }
    uint32_t hashTimeAccu = 0, bfTimeAccu = 0;
    chrono::time_point<chrono::high_resolution_clock> t1, t2, t3;
    // # pragma omp parallel for 
    for (size_t i = 0; i < sequences.size(); ++i) {
        uint64_t hashes[config.k];
        uint32_t threadId = 0; // omp_get_thread_num();
        hasher[threadId]->setSequence(sequences[i]);
        while (hasher[threadId]->hasNext()) {
            t1 = chrono::high_resolution_clock::now();
            hasher[threadId]->hash(hashes);
            t2 = chrono::high_resolution_clock::now();
            bf.insert(hashes);
            t3 = chrono::high_resolution_clock::now();
            hashTimeAccu += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
            bfTimeAccu += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
        }
    }
    cout << "Hash Time: " << hashTimeAccu << "; Insert Time: " << bfTimeAccu << endl;
    cout << "BF Packing: " << bf.count() << '/' << bf.size << endl;
    cout << "# of items inserted: " << sequences.size() << endl;

    // hashTimeAccu = 0, bfTimeAccu = 0;
    uint32_t fpCount = 0;
    // # pragma omp parallel for 
    for (size_t i = 0; i < sequences.size(); ++i) {
        uint32_t threadId = 0; // omp_get_thread_num();
        uint64_t hashes[config.k];
        sequences[i] = poison(i, sequences[i], config.poison);
        hasher[threadId]->setSequence(sequences[i]);
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
    cout << "Hash Time: " << hashTimeAccu << "; Query Time: " << bfTimeAccu << endl;
    float fpRate = static_cast<float>(fpCount) / sequences.size();
    cout << "False Positive Rate: " << fpRate << endl;
    fprintf(stderr, "%d / %d / %.2e\n", hashTimeAccu, bfTimeAccu, fpRate);

    bf.release();
    return 0;
}
