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

class FuzzyHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    const uint32_t minHashRange;
public:
    FuzzyHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange), minHashRange(range / universalHashRange) {
        if (range % universalHashRange != 0) {
            cerr << "Hash range=" << range << " is not a multiple of " << universalHashRange << '.' << endl;
        }
    }

    void hash(uint32_t * out) override {
        for (uint32_t i = 0; i < num_hashes; ++i) {
            uint32_t hashValue, minValue = UINT32_MAX;
            for (uint32_t j = 0; j <= 31 - kMer; ++j) {
                MurmurHash3_x86_32(sequence + pos + j, kMer, seed + i * 31 + j, &hashValue);
                minValue = min(hashValue, minValue);
            }
            MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, universalHashRange);
            out[i] += (minValue % minHashRange) * universalHashRange;
        }
        pos += 1;
    }
};

class EfficientFuzzyHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    const uint32_t minHashRange;
    uint32_t *tree;
    uint8_t *indices;
public:
    EfficientFuzzyHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange), minHashRange(range / universalHashRange), tree(new uint32_t[63 * 4]), indices(new uint8_t[31 * 4]) {
        if (range % universalHashRange != 0) {
            cerr << "Hash range=" << range << " is not a multiple of " << universalHashRange << '.' << endl;
        }
        fill(tree, tree + 63 * 4, UINT32_MAX);
        for (uint32_t i = 0; i < 31 * 4; ++i) {
            indices[i] = i % (31 * 4);
        }
    }

    void hash(uint32_t * out) override {
        if (tree[0] == UINT32_MAX) {
            // initial tree construction
            for (uint32_t i = 0; i < num_hashes; ++i) {
                uint32_t hashValue, minValue = UINT32_MAX;
                for (uint32_t j = 0; j <= 31 - kMer; ++j) {
                    MurmurHash3_x86_32(sequence + pos + j, kMer, seed + i * 31 + j, &hashValue);
                    minValue = min(hashValue, minValue);
                }
                MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, universalHashRange);
                out[i] += (minValue % minHashRange) * universalHashRange;
            }
            pos += 1;
        }
    }
};

#endif
