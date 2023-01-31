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
#include <filesystem>

using namespace std;

int main(int argc, char** argv) {
    const Config config = getConfigs(argv[1], argc - 2, argv + 2);
    config.print();

    if (!filesystem::exists(config.experimentDir)) {
        filesystem::create_directory(config.experimentDir);
    }

    ofstream results(config.experimentDir + '/' + config.resultFile, ofstream::trunc);
    config.print(results);

    string filellistname = config.fastqFileName;
    bool insert = !config.queryOnly;
    bool query = true;
    // Create and open a log file
    // ofstream InsertLogs("./results/insertLogs.txt");
    Hasher *hasher[config.numThreads]; // each thread gets its own hasher
    for (uint32_t i = 0; i < config.numThreads; ++i) {
        switch(config.hashType) {
            case Config::MURMUR_HASH:
                hasher[i] = static_cast<Hasher*>(new MurmurHasher(UINT64_MAX, config.k, config.seed));
                break;
            case Config::BRUTE_FORCE_FUZZY_HASH:
                hasher[i] = static_cast<Hasher*>(new FuzzyHasher(UINT64_MAX, config.k, config.kMer, config.universalHashRange, config.seed));
                break;
            case Config::FUZZY_HASH:
                hasher[i] = static_cast<Hasher*>(new EfficientFuzzyHasher(UINT64_MAX, config.k, config.kMer, config.universalHashRange, config.seed));
                break;
            case Config::FUZZY_HASH_EXP:
                hasher[i] = static_cast<Hasher*>(new EfficientFuzzyHasherEXP(UINT64_MAX, config.k, config.kMer, config.universalHashRange, config.seed));
                break;
            default:
                cout << "HASH TYPE not recognized" << endl;
                return 0;
        }
    }

    BloomFilter** ArBF_array;
    uint32_t numFiles = config.firstNFiles;
    ArBF_array = new BloomFilter*[numFiles]; //array of pointers
    vector<pair<uint32_t, uint32_t>> fileStats(numFiles, make_pair(0, 0));
    vector<uint64_t> ranges;
    uint64_t range;

    vector<string> fileNames(numFiles);

    if (insert) {
        uint32_t n = 0, hashTimeAccu = 0, bfTimeAccu = 0, counter = 0;
        chrono::time_point<chrono::high_resolution_clock> t1, t2, t3;
        ifstream filellist(filellistname);
        while (n < numFiles && getline(filellist, fileNames[n])) {
            vector<string> sequences = getFastqData("/scratch/gg29/CKBF/data/fastqFiles/" + fileNames[n]);
            uint64_t numInserts = 0;
            for (size_t s = 0; s < sequences.size(); ++s) {
                numInserts += sequences[s].size() - 30;
            }
            range = static_cast<uint64_t>(-config.rangefactor * numInserts * config.k / log(1 - pow(config.fpr, 1.0 / config.k)));
            // range = static_cast<uint64_t>(config.rangefactor * config.k * numInserts);
            // range = static_cast<uint64_t>(config.rangefactor) * static_cast<uint64_t>(config.k) * static_cast<uint64_t>(sequences.size());
            ArBF_array[n] = new BloomFilter(range, config.k, config.disk, false, config.experimentDir + '/' + fileNames[n].substr(0, fileNames[n].length() - 6) + ".dat");
            uint32_t hashTime = 0, bfTime = 0;
            // # pragma omp parallel for 
            for (size_t i = 0; i < sequences.size(); ++i) {
                uint64_t hashes[config.k];
                uint32_t threadId = 0; // omp_get_thread_num();
                hasher[threadId]->setSequence(sequences[i]);
                while (hasher[threadId]->hasNext()) {
                    t1 = chrono::high_resolution_clock::now();
                    hasher[threadId]->hash(hashes);
                    for (uint8_t k = 0; k < config.k; ++k){
                        hashes[k] = hashes[k] % range;
                    }
                    t2 = chrono::high_resolution_clock::now();
                    ArBF_array[n]->insert(hashes);
                    t3 = chrono::high_resolution_clock::now();
                    hashTime += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
                    bfTime += chrono::duration_cast<chrono::microseconds>(t3 - t2).count();
                    counter += 1;
                }
                if (100 * (i-1) / sequences.size() != 100 * i / sequences.size()) {
                    cout << '#' << flush;
                }
            }
            cout << "\nBF #" << n << ' ' << fileNames[n] << " Packing: " << ArBF_array[n]->count() << '/' << ArBF_array[n]->size
                << "; # of inserts: " << counter << "; Hash Time: " << hashTime
                << "; BF Insert Time: " << bfTime << endl;
            results << "BF #" << n << ' ' << fileNames[n] << " Packing: " << ArBF_array[n]->count() << '/' << ArBF_array[n]->size
                << "; # of inserts: " << counter << "; Hash Time: " << hashTime
                << "; BF Insert Time: " << bfTime << endl;

            // ArBF_array[n]->release();
            n++;
            hashTimeAccu += hashTime;
            bfTimeAccu += bfTime;
            ranges.push_back(range);
        }
        if (n < numFiles) {
            numFiles = n;
        }
        cout << "Indexed " << numFiles << " .fastq files." << endl;
        cout << "Total Hash Time: " << hashTimeAccu << "; Total BF Insert Time: " << bfTimeAccu << endl;
        results << "Total Hash Time: " << hashTimeAccu << "; Total BF Insert Time: " << bfTimeAccu << endl;
        write_vector(config.experimentDir + "/ranges.txt", ranges);
    }

    if (query) {
        uint32_t n = 0;
        string queryFile = config.queryFileName;
        vector<string> queries;
        getQueryForArBF(queryFile, queries);
        uint32_t nq = queries.size();

        uint32_t hashTimeAccu = 0, bfTimeAccu = 0;
        chrono::time_point<chrono::high_resolution_clock> t1, t2;

        // ranges = read_vector(config.experimentDir + "/ranges.txt");
        // ifstream filellist(filellistname);
        // while (n < numFiles && getline(filellist, fileNames[n])) {
        //     range = ranges[n];
        //     ArBF_array[n] = new BloomFilter(range, config.k, config.disk, true, config.experimentDir + '/' + fileNames[n].substr(0,fileNames[n].length() - 6) + ".dat");
        //     n++;
        // }
        uint32_t nPositives = 0;

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

            // results << queries[i];

            for (size_t n = 0; n < numFiles; ++n){
                bool positive = true;
                t1 = chrono::high_resolution_clock::now();
                for (uint32_t l = 0; l < nHashes; ++l) {
                    positive = positive && ArBF_array[n]->test(hashes + l * config.k, ranges[n]);
                }
                t2 = chrono::high_resolution_clock::now();
                uint32_t duration = chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
                fileStats[n].first += duration;
                bfTimeAccu += duration;
                if (positive) {
                    // results << ',' << n;
                    fileStats[n].second += 1;
                    nPositives += 1;
                }
            }
            // results << endl;
        }

        for (uint32_t i = 0; i < numFiles; ++i) {
            results << "BF #" << i << ' ' << fileNames[i] << " Lookup Time: " << fileStats[i].first << "; # of positives: " << fileStats[i].second << endl;
            // ArBF_array[i]->release();
        }

        cout << "Query Hash Time: " << hashTimeAccu << "; BF Lookup Time: " << bfTimeAccu << "; # of positives: " << nPositives << '/' << nq * numFiles << endl;
        results << "Query Hash Time: " << hashTimeAccu << "; BF Lookup Time: " << bfTimeAccu << "; # of positives: " << nPositives << '/' << nq * numFiles << endl;
    }
    return 0;
}

//todo: 
// 1) choice of ABF experiment in RAM, use disk flag (if 0) to save and load BF arrays n RAM 
// 2) stop when first kmer gives False, or calculate the %  kmer content
// 3) get query,fileID queryfile