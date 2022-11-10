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
    const uint32_t cacheLineSize = 512;
    const uint32_t minHashRange;
    const uint32_t kMer;
public:
    FuzzyHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), minHashRange(range / cacheLineSize), kMer(kMer) {
        if (range % cacheLineSize != 0) {
            cerr << "Hash range=" << range << " is not a multiple of " << cacheLineSize << '.' << endl;
        }
    }

    void hash(uint32_t * out) override {
        for (uint32_t i = 0; i < num_hashes; ++i) {
            uint32_t hashValue, minValue = UINT32_MAX;
            for (uint32_t j = 0; j <= 31 - kMer; ++j) {
                MurmurHash3_x86_32(sequence + pos + j, kMer, seed + i * 31 + j, &hashValue);
                minValue = min(hashValue, minValue);
            }
            MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, cacheLineSize);
            out[i] += (minValue % minHashRange) * cacheLineSize;
        }
        pos += 1;
    }
};

#endif
