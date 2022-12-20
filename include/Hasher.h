#ifndef _HASHER_H_
#define _HASHER_H_

#include <iostream>
#include <string>
#include "MurmurHash3.h"

using namespace std;

class Hasher {
protected:
    const char *sequence;
    size_t sequence_len;
    const uint32_t range;
    const uint32_t num_hashes;
    const uint32_t seed;
    uint32_t pos;
public:
    Hasher(uint32_t range, uint32_t num_hashes, uint32_t seed=0)
        : range(range), num_hashes(num_hashes), seed(seed), pos(0) {}

    inline bool hasNext() {
        return pos <= sequence_len - 31;
    }

    virtual void hash(uint32_t * out) = 0;

    void setSequence(const string & _sequence) {
        sequence = _sequence.data();
        sequence_len = _sequence.size();
        pos = 0;
    }
};

class MurmurHasher : public Hasher {
public:
    MurmurHasher(uint32_t range, uint32_t num_hashes, uint32_t seed=0)
        : Hasher(range, num_hashes, seed) {}

    void hash(uint32_t * out) override {
        for (uint32_t j = 0; j < num_hashes; ++j) {
            MurmurHash3_x86_32(sequence + pos, 31, seed + j, out + j, range);
        }
        pos += 1;
    }
};


class LSHBBFHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    const uint32_t minHashRange;
public:
    LSHBBFHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange), minHashRange(range / universalHashRange) {
        if (range % universalHashRange != 0) {
            cerr << "Hash range=" << range << " is not a multiple of " << universalHashRange << '.' << endl;
        }
    }

    void hash(uint32_t * out) override {
        for (uint32_t i = 0; i < num_hashes; ++i) {
            uint32_t hashValue, minValue = UINT32_MAX;
            for (uint32_t j = 0; j <= 31 - kMer; ++j) {
                MurmurHash3_x86_32(sequence + pos + j, kMer, seed + i, &hashValue);
                minValue = min(hashValue, minValue);
            }
            MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, universalHashRange);
            out[i] += (minValue % minHashRange) * universalHashRange;
        }
        pos += 1;
    }
};

class FuzzyHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    const uint32_t minHashRange;
public:
    FuzzyHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange), minHashRange(range) {
    }

    void hash(uint32_t * out) override {
        uint32_t hashValue, minValue = UINT32_MAX;
        for (uint32_t i = 0; i < num_hashes; ++i) {
            minValue = UINT32_MAX;
            for (uint32_t j = 0; j <= 31 - kMer; ++j) {
                MurmurHash3_x86_32(sequence + pos + j, kMer, seed + i, &hashValue);
                minValue = min(hashValue, minValue);
            }
            MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, universalHashRange);
            // offset with min hash. overall uh + mh
            out[i] += minValue % (minHashRange - universalHashRange); 
        }
        pos += 1;
    }
};

class EfficientFuzzyHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    const uint32_t minHashRange;
    uint32_t *trees;
    uint32_t index_to_pop;
public:
    EfficientFuzzyHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange), minHashRange(range), trees(new uint32_t[((32 - kMer) * 2 - 1) * num_hashes]), index_to_pop(0) {
        fill(trees, trees + ((32 - kMer) * 2 - 1) * num_hashes, UINT32_MAX);
    }

    void buildTrees() {
        for (uint32_t i = 0; i < num_hashes; ++i) {
            for (int32_t j = kMer - 2; j >= 0; --j) {
                trees[((32 - kMer) * 2 - 1) * i + j] = min(
                    trees[((32 - kMer) * 2 - 1) * i + j * 2 + 1],
                    trees[((32 - kMer) * 2 - 1) * i + j * 2 + 2]
                );
            }
        }
    }

    void hash(uint32_t * out) override {
        if (pos == 0) {
            // initial tree construction
            uint32_t hashValue, minValue;
            for (uint32_t i = 0; i < num_hashes; ++i) {
                minValue = UINT32_MAX;
                for (uint32_t j = 0; j <= 31 - kMer; ++j) {
                    MurmurHash3_x86_32(sequence + pos + j, kMer, seed + i, &hashValue);
                    minValue = min(hashValue, minValue);
                    trees[((32 - kMer) * 2 - 1) * i + (kMer - 1 + j)] = hashValue;
                }
                MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, universalHashRange);
                out[i] += minValue % (minHashRange - universalHashRange);
            }
            buildTrees();
        } else {
            uint32_t hashValue, new_index;
            for (uint32_t i = 0; i < num_hashes; ++i) {
                new_index = kMer - 1 + index_to_pop;
                MurmurHash3_x86_32(sequence + pos + 31 - kMer, kMer, seed + i, &hashValue);
                trees[((32 - kMer) * 2 - 1) * i + new_index] = hashValue;
                while (new_index > 0) {
                    new_index = (new_index - 1) / 2;
                    trees[((32 - kMer) * 2 - 1) * i + new_index] = min(
                        trees[((32 - kMer) * 2 - 1) * i + new_index * 2 + 1],
                        trees[((32 - kMer) * 2 - 1) * i + new_index * 2 + 2]
                    );
                }
                MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, universalHashRange);
                out[i] += trees[((32 - kMer) * 2 - 1) * i] % (minHashRange - universalHashRange);
            }
            index_to_pop = (index_to_pop + 1) % (32 - kMer);
        }
        pos += 1;
    }
};

#endif
