#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <iostream>
#include <string>
#include <fstream>

using namespace std;

struct Config {
    static const uint32_t MURMUR_HASH = 0;
    static const uint32_t FUZZY_HASH = 1;
    static const uint32_t BRUTE_FORCE_FUZZY_HASH = 2;

    string fastqFileName;
    string queryFileName;
    uint64_t range;
    uint32_t k; // number of hashes
    uint32_t numThreads;
    uint32_t universalHashRange = 0;
    uint32_t seed;
    uint32_t hashType;
    uint32_t kMer = 0;
    uint32_t poison = 0;
    bool disk = false;

    void print() const {
        cout << "configs: fastqFileName=" << fastqFileName << "; queryFileName=" << queryFileName << "; range=" << range << "; k=" << k << "; numThreads=" << numThreads << "; universalHashRange=" << universalHashRange << "; kMer=" << kMer << "; seed=" << seed << "; hashType=" << (hashType ? "fuzzy" : "murmur") << "; poison=" << poison << "; disk=" << (disk ? "yes" : "no") << endl;
    }
};

Config getConfigs(string configFileName) {
    ifstream configFile(configFileName);
    Config config;

    string line;
    while(getline(configFile, line)) {
        size_t equalPos = line.find('=');
        string key = line.substr(0, equalPos);
        string value = line.substr(equalPos + 1);
        if (key == "fastq_file_name") {
            config.fastqFileName = value;
        } else if (key == "query_file_name") {
            config.queryFileName = value;
        } else if (key == "range") {
            config.range = stoul(value);
        } else if (key == "k") {
            config.k = stoul(value);
        } else if (key == "num_threads") {
            config.numThreads = stoul(value);
        } else if (key == "universal_hash_range") {
            config.universalHashRange = stoul(value);
        } else if (key == "kmer") {
            config.kMer = stoul(value);
        } else if (key == "seed") {
            config.seed = stoul(value);
        } else if (key == "hash") {
            if (value == "murmur") {
                config.hashType = Config::MURMUR_HASH;
            } else if (value == "fuzzy") {
                config.hashType = Config::FUZZY_HASH;
            } else if (value == "brute_force_fuzzy") {
                config.hashType = Config::BRUTE_FORCE_FUZZY_HASH;
            } else {
                cerr << value << " is an unrecognized hash type" << endl;
            }
        } else if (key == "poison") {
            config.poison = stoul(value);
        } else if (key == "disk") {
            config.disk = stoul(value) > 0;
        } else {
            cerr << key << " is an unrecognized config key" << endl;
        }
    }
    return config;
}

#endif
