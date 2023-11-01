#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <iostream>
#include <string>
#include <fstream>

using namespace std;

struct Config {
    string hashtypeStr[5]  = {
      "MURMUR_HASH",
      "FUZZY_HASH",
      "BRUTE_FORCE_FUZZY_HASH",
      "FUZZY_HASH_EXP",
      "MAX_HASH_TYPE"
    };

    enum HASHTYPE{
      MURMUR_HASH = 0,
      FUZZY_HASH = 1,
      BRUTE_FORCE_FUZZY_HASH = 2,
      FUZZY_HASH_EXP= 3,
      MAX_HASH_TYPE = 4
    };

    string dataPath;
    string fastqFileName;
    string queryFileName;
    uint64_t range;
    uint32_t k; // number of hashes
    uint32_t numThreads;
    uint32_t universalHashRange = 0;
    uint32_t seed;
    HASHTYPE hashType = MAX_HASH_TYPE;
    uint32_t kMer = 0;
    uint32_t poison = 0;
    uint32_t firstNFiles = 200;
    bool queryOnly = false;
    bool disk = false;
    float rangefactor=1;

    void print() const {
        cout << "configs: fastqFileName=" << fastqFileName << "; queryFileName=" << queryFileName 
        << "; range=" << range << "; k=" << k << "; numThreads=" << numThreads 
        << "; universalHashRange=" << universalHashRange << "; kMer=" << kMer << "; seed=" 
        << seed << "; hashType=" << hashtypeStr[hashType] << "; poison=" << poison 
        << "; disk=" << (disk ? "yes" : "no") << "; rangefactor=" << rangefactor
        << "; data_path=" << dataPath
        << "; first_n_files=" << firstNFiles << "; query_only=" << (queryOnly ? "yes" : "no") << endl;
    }
};

Config getConfigs(string configFileName, int numOverrides = 0, char ** overrideConfigs = nullptr) {
    ifstream configFile(configFileName);
    Config config;

    string line;
    while (true) {
        if (!getline(configFile, line)) {
            if (numOverrides > 0) {
                line = string(overrideConfigs[--numOverrides]);
            } else {
                break;
            }
        }

        size_t equalPos = line.find('=');
        string key = line.substr(0, equalPos);
        string value = line.substr(equalPos + 1);
        if (key == "fastq_file_name") {
            config.fastqFileName = value;
        } else if (key == "query_file_name") {
            config.queryFileName = value;
        } else if (key == "data_path") {
            config.dataPath = value;
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
            } else if (value == "fuzzyexp") {
                config.hashType = Config::FUZZY_HASH_EXP;
            } else if (value == "fuzzy") {
                config.hashType = Config::FUZZY_HASH;
            } else if (value == "brute_force_fuzzy") {
                config.hashType = Config::BRUTE_FORCE_FUZZY_HASH;
            } else {
                cerr << value << " is an unrecognized hash type" << endl;
                config.hashType = Config::MAX_HASH_TYPE;
            }
        } else if (key == "poison") {
            config.poison = stoul(value);
        } else if (key == "disk") {
            config.disk = stoul(value) > 0;
        } else if (key == "rangefactor") {
            config.rangefactor = stoul(value);
        } else if (key == "first_n_files") {
            config.firstNFiles = stoul(value);
        } else if (key == "query_only") {
            config.queryOnly = stoul(value) > 0;
        } else {
            cerr << key << " is an unrecognized config key" << endl;
        }
    }
    return config;
}

#endif
