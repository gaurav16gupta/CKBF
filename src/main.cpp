#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "BloomFilter.h"
#include "utils.h"
#include "config.h"
#include "Hasher.h"

using namespace std;

int main() {
    Config config = getConfigs("../configs/default.cfg");
    config.print();

    vector<string> sequences = getFastqData("../data/" + config.fastqFileName + ".fastq");
    BloomFilter bf(config.range, config.k);
    Hasher *hasher = config.hashType == Config::MURMUR_HASH
        ? static_cast<Hasher*>(new MurmurHasher(config.range, config.k, config.seed))
        : static_cast<Hasher*>(new FuzzyHasher(config.range, config.k, 3, config.seed));

    for (size_t i = 0; i < sequences.size(); ++i) {
        uint32_t hashes[config.k];
        hasher->setSequence(sequences[i]);
        while (hasher->hasNext()) {
            hasher->hash(hashes);
            bf.insert(hashes);
        }
    }

    uint32_t fp = 0;
    for (size_t i = 0; i < sequences.size(); ++i) {
        sequences[i][10] = 'B';
        uint32_t hashes[config.k];
        hasher->setSequence(sequences[i]);
        bool pos = true;
        while (hasher->hasNext()) {
            hasher->hash(hashes);
            pos &= bf.test(hashes);
        }
        fp += pos;
    }

    cout << fp << endl;

    return 0;
}