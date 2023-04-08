#ifndef _HASHER_H_
#define _HASHER_H_

#include <iostream>
#include <string>
#include "MurmurHash3.h"
#include <assert.h>
#include <cmath>
using namespace std;

class Hasher {
protected:
    const char *sequence;
    size_t sequence_len;
    const uint64_t range;
    const uint32_t num_hashes;
    const uint32_t seed;
    uint32_t pos;
    uint64_t RANGEMASK;
public:
    Hasher(uint64_t range, uint32_t num_hashes, uint32_t seed=0)
        : range(range), num_hashes(num_hashes), seed(seed), pos(0) {

        uint64_t Rbits = (int)std::log2((float) range);
        RANGEMASK = (1LL<<Rbits) -1;
    }

    inline bool hasNext() {
        return pos <= sequence_len - 31;
    }

    uint32_t nHashesInSequence() {
        return sequence_len - 30;
    }

    virtual void hash(uint64_t * out) = 0;

    void setSequence(const string & _sequence) {
        sequence = _sequence.data();
        sequence_len = _sequence.size();
        pos = 0;
    }
};

class MurmurHasher : public Hasher {
public:
    MurmurHasher(uint64_t range, uint32_t num_hashes, uint32_t seed=0)
        : Hasher(range, num_hashes, seed) {}

    void hash(uint64_t * out) override {
        for (uint32_t j = 0; j < num_hashes; ++j) {
            MurmurHash3_x64_64(sequence + pos, 31, seed + j, out + j);
            out[j] = out[j]  & RANGEMASK;
        }
        pos += 1;
    }
};


// class LSHBBFHasher : public Hasher {
// protected:
//     const uint32_t kMer;
//     const uint32_t universalHashRange;
//     const uint32_t minHashRange;
// public:
//     LSHBBFHasher(uint32_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
//         : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange), minHashRange(range / universalHashRange) {
//         if (range % universalHashRange != 0) {
//             cerr << "Hash range=" << range << " is not a multiple of " << universalHashRange << '.' << endl;
//         }
//     }

//     void hash(uint32_t * out) override {
//         for (uint32_t i = 0; i < num_hashes; ++i) {
//             uint32_t hashValue, minValue = UINT32_MAX;
//             for (uint32_t j = 0; j <= 31 - kMer; ++j) {
//                 MurmurHash3_x86_32(sequence + pos + j, kMer, seed + i, &hashValue);
//                 minValue = min(hashValue, minValue);
//             }
//             MurmurHash3_x86_32(sequence + pos, 31, seed + i, out + i, universalHashRange);
//             out[i] += (minValue % minHashRange) * universalHashRange;
//         }
//         pos += 1;
//     }
// };

class FuzzyHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    uint64_t URANGEMASK;
public:
    FuzzyHasher(uint64_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange) {
        uint64_t Lbits = (int)std::log2((float) universalHashRange);
        URANGEMASK = (1LL<<Lbits) -1;
    }

    void hash(uint64_t * out) override {
        uint64_t hashValue, minValue;
        for (uint32_t i = 0; i < num_hashes; ++i) {
            minValue = UINT64_MAX;
            for (uint32_t j = 0; j <= 31 - kMer; ++j) {
                MurmurHash3_x64_64(sequence + pos + j, kMer, seed + i, &hashValue);
                minValue = min(hashValue, minValue);
            }
            MurmurHash3_x64_64(sequence + pos, 31, seed + i, out + i);
            out[i] = out[i]  & URANGEMASK;
            // offset with min hash. overall uh + mh
            out[i] += minValue & RANGEMASK; 
        }
        pos += 1;
    }
};

static void print_trees(uint64_t *trees, uint32_t len, uint32_t num) {
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

class EfficientFuzzyHasherEXP : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    uint32_t tree_length;
    uint64_t *trees;
    uint32_t index_to_pop;
    uint64_t* oneperm_hash_helper_arr = new uint64_t [33];
    int32_t Lbits;
    uint64_t URANGEMASK;
public:
    EfficientFuzzyHasherEXP(uint64_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange),
          tree_length(next_pow_of_2(32 - kMer) * 2 - 1), trees(new uint64_t[tree_length * num_hashes]), index_to_pop(0) {
        fill(trees, trees + tree_length * num_hashes, UINT64_MAX);
        Lbits = (int)std::log2((float) universalHashRange);
        URANGEMASK = (1LL<<Lbits) -1;
        std::cout << "Lbits: " << Lbits << "urange: " << universalHashRange << std::flush << std::endl;
        assert(num_hashes*Lbits <= 64);
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

    void hash(uint64_t * out) override {
        uint64_t hashValue, minValue, uhashvalue;
        uint32_t i, j, new_index;
        uint32_t idx1;
        //uint32_t seed1;
        uint64_t one_perm_shift = UINT64_MAX / num_hashes;
        MurmurHash3_x64_64(sequence + pos, 31, seed, &uhashvalue);
        if (pos == 0) {
            // initial tree construction
            for (j=0; j<=32 - kMer; ++j) {
                MurmurHash3_x64_64(sequence + pos + j, kMer, seed, &oneperm_hash_helper_arr[j]);
            }

            for (i = 0; i < num_hashes; ++i) {
                minValue = UINT64_MAX;
                idx1 = tree_length * i + tree_length / 2;
                //seed1 = seed + i;
                for (j = 0; j <= 31 - kMer; ++j) {
                    hashValue = oneperm_hash_helper_arr[j] - i*one_perm_shift;
                    minValue = min(hashValue, minValue);
                    trees[idx1 + j] = hashValue;
                }
                //MurmurHash3_x64_64(sequence + pos, 31, seed1, out + i, universalHashRange);
                out[i] = (uhashvalue >> (i*Lbits)) & URANGEMASK;
                out[i] += minValue & RANGEMASK;
            }
            buildTrees();
            index_to_pop = 0;
        } else {
            MurmurHash3_x64_64(sequence + pos + 31 - kMer, kMer, seed, &oneperm_hash_helper_arr[0]);
            for (i = 0; i < num_hashes; ++i) {
                hashValue = (oneperm_hash_helper_arr[0] -i*one_perm_shift);
                new_index = tree_length / 2 + index_to_pop;
                //seed1 = seed + i;
                idx1 = tree_length * i;
                trees[idx1 + new_index] = hashValue;
                while (new_index > 0) {
                    new_index = (new_index - 1) / 2;
                    trees[idx1 + new_index] = min(
                        trees[idx1 + new_index * 2 + 1],
                        trees[idx1 + new_index * 2 + 2]
                    );
                }
                //MurmurHash3_x64_64(sequence + pos, 31, seed1, out + i, universalHashRange);
                out[i] = (uhashvalue >> (i*Lbits)) & URANGEMASK;
                out[i] += trees[idx1] & RANGEMASK;
            }
            index_to_pop = (index_to_pop + 1) % (32 - kMer);
        }
        pos += 1;
    }
};

class EfficientFuzzyHasher : public Hasher {
protected:
    const uint32_t kMer;
    const uint32_t universalHashRange;
    uint32_t tree_length;
    uint64_t *trees;
    uint32_t index_to_pop;
    uint64_t URANGEMASK;
public:
    EfficientFuzzyHasher(uint64_t range, uint32_t num_hashes, uint32_t kMer, uint32_t universalHashRange, uint32_t seed=0)
        : Hasher(range, num_hashes, seed), kMer(kMer), universalHashRange(universalHashRange),
          tree_length(next_pow_of_2(32 - kMer) * 2 - 1), trees(new uint64_t[tree_length * num_hashes]), index_to_pop(0) {
        fill(trees, trees + tree_length * num_hashes, UINT64_MAX);
        uint64_t Lbits = (int)std::log2((float) universalHashRange);
        URANGEMASK = (1LL<<Lbits) -1;
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

    void hash(uint64_t * out) override {
        uint64_t hashValue, minValue;
        uint32_t i, j, new_index;
        uint32_t idx1, seed1;
        if (pos == 0) {
            // initial tree construction
            for (i = 0; i < num_hashes; ++i) {
                minValue = UINT64_MAX;
                idx1 = tree_length * i + tree_length / 2;
                seed1 = seed + i;
                for (j = 0; j <= 31 - kMer; ++j) {
                    MurmurHash3_x64_64(sequence + pos + j, kMer, seed1, &hashValue);
                    minValue = min(hashValue, minValue);
                    trees[idx1 + j] = hashValue;
                }
                MurmurHash3_x64_64(sequence + pos, 31, seed1, out + i);
                out[i] = out[i] & URANGEMASK;
                out[i] += minValue & RANGEMASK;
            }
            buildTrees();
            index_to_pop = 0;
        } else {
            for (i = 0; i < num_hashes; ++i) {
                new_index = tree_length / 2 + index_to_pop;
                seed1 = seed + i;
                idx1 = tree_length * i;
                MurmurHash3_x64_64(sequence + pos + 31 - kMer, kMer, seed1, &hashValue);
                trees[idx1 + new_index] = hashValue;
                while (new_index > 0) {
                    new_index = (new_index - 1) / 2;
                    trees[idx1 + new_index] = min(
                        trees[idx1 + new_index * 2 + 1],
                        trees[idx1 + new_index * 2 + 2]
                    );
                }
                MurmurHash3_x64_64(sequence + pos, 31, seed1, out + i);
                out[i] = out[i] & URANGEMASK;
                out[i] += trees[idx1] & RANGEMASK;
            }
            index_to_pop = (index_to_pop + 1) % (32 - kMer);
        }
        pos += 1;
    }

};

#endif
