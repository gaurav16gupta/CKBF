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

static void print_trees(uint32_t *trees, uint32_t len, uint32_t num) {
    for (uint32_t j = 0; j < num; ++j) {
        for (uint32_t i = 0; i < len; ++i) {
            cout << trees[len * j + i] << ' ';
            if (!((i + 1) & (i + 2))) {
                cout << endl;
            }
        }
        cout << endl;
    }
    cout << endl;
}

static uint32_t next_pow_of_2(uint32_t x) {
    uint32_t power = 1;
    while (power < x) {
        power *= 2;
    }
    return power;
}

class EfficientFuzzyHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    const uint32_t minHashRange;
    uint32_t tree_length;
    uint32_t *trees;
    uint32_t index_to_pop;
public:
    EfficientFuzzyHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange), minHashRange(range),
          tree_length(next_pow_of_2(32 - kMer) * 2 - 1), trees(new uint32_t[tree_length * num_hashes]), index_to_pop(0) {
        fill(trees, trees + tree_length * num_hashes, UINT32_MAX);
    }

    void buildTrees() {
        uint32_t idx1;
        for (uint32_t i = 0; i < num_hashes; ++i) {
            idx1 = tree_length * i;
            for (int32_t j = tree_length / 2 - 1; j >= 0; --j) {
                trees[idx1 + j] = min(
                    trees[idx1 + j * 2 + 1],
                    trees[idx1 + j * 2 + 2]
                );
            }
        }
    }

    void hash(uint32_t * out) override {
        uint32_t hashValue, minValue, i,j, new_index;
        uint32_t idx1, seed1;
        if (pos == 0) {
            // initial tree construction
            for (i = 0; i < num_hashes; ++i) {
                minValue = UINT32_MAX;
                idx1 = tree_length * i + tree_length / 2;
                seed1 = seed + i;
                for (j = 0; j <= 31 - kMer; ++j) {
                    MurmurHash3_x86_32(sequence + pos + j, kMer, seed1, &hashValue);
                    minValue = min(hashValue, minValue);
                    trees[idx1 + j] = hashValue;
                }
                MurmurHash3_x86_32(sequence + pos, 31, seed1, out + i, universalHashRange);
                out[i] += minValue % (minHashRange - universalHashRange);
            }
            buildTrees();
            index_to_pop = 0;
        } else {
            for (i = 0; i < num_hashes; ++i) {
                new_index = tree_length / 2 + index_to_pop;
                seed1 = seed + i;
                idx1 = tree_length * i;
                MurmurHash3_x86_32(sequence + pos + 31 - kMer, kMer, seed1, &hashValue);
                trees[idx1 + new_index] = hashValue;
                while (new_index > 0) {
                    new_index = (new_index - 1) / 2;
                    trees[idx1 + new_index] = min(
                        trees[idx1 + new_index * 2 + 1],
                        trees[idx1 + new_index * 2 + 2]
                    );
                }
                MurmurHash3_x86_32(sequence + pos, 31, seed1, out + i, universalHashRange);
                out[i] += trees[idx1] % (minHashRange - universalHashRange);
            }
            index_to_pop = (index_to_pop + 1) % (32 - kMer);
        }
        pos += 1;
    }

};

class EfficientOnePermFuzzyHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    const uint32_t minHashRange;
    uint32_t *trees;
    uint32_t index_to_pop;
public:
    EfficientOnePermFuzzyHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange), minHashRange(range), trees(new uint32_t[((32 - kMer) * 2 - 1) * num_hashes]), index_to_pop(0) {
        fill(trees, trees + ((32 - kMer) * 2 - 1) * num_hashes, UINT32_MAX);
    }

    void buildTrees() {
        uint32_t idx1;
        for (uint32_t i = 0; i < num_hashes; ++i) {
            idx1 = ((32 - kMer) * 2 - 1) * i ;
            for (int32_t j = kMer - 2; j >= 0; --j) {
                trees[idx1 + j] = min(
                    trees[idx1 + j * 2 + 1],
                    trees[idx1 + j * 2 + 2]
                );
            }
        }
    }

    void hash(uint32_t * out) {
        uint32_t hashValue, i,j, new_index, idx1, this_hashvalue;
        uint32_t stride = UINT32_MAX / num_hashes;
        if (pos == 0) {
            idx1 = ((32 - kMer) * 2 - 1);
            for(j=0; j<=31-kMer; ++j) {
                MurmurHash3_x86_32(sequence + pos + j, kMer, seed, &hashValue);
                for (i=0; i<=num_hashes; ++i) {
                    this_hashvalue = hashValue - i*stride; // this is shifting and taking mod with 2^32
                    trees[idx1 * i +  (kMer - 1) +  j] = this_hashvalue;
                }
            }
            buildTrees();
            for(i=0;i<=num_hashes;++i) {
                MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, universalHashRange);
                out[i] += trees[idx1*i] % (minHashRange - universalHashRange);
            }
        } else {
            MurmurHash3_x86_32(sequence + pos + 31 - kMer, kMer, seed, &hashValue);
            for (i = 0; i < num_hashes; ++i) {
                this_hashvalue = hashValue - i*stride; // this is shifting and taking mod with 2^32
                new_index = kMer - 1 + index_to_pop;
                idx1 = ((32 - kMer) * 2 - 1) * i;
                trees[idx1 + new_index] = this_hashvalue;
                while (new_index > 0) {
                    new_index = (new_index - 1) / 2;
                    trees[idx1 + new_index] = min(
                        trees[idx1 + new_index * 2 + 1],
                        trees[idx1 + new_index * 2 + 2]
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
