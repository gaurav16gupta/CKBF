#ifndef _MYBLOOM_
#define _MYBLOOM_
#include <vector>
#include <bitset>
#include "bitArray.h"


class BloomFilter{
    public:
        // BloomFilter(int capacity, float FPR, int k);
        BloomFilter(int sz, int k);
        void insert(std::vector<uint> a);
        bool test(std::vector<uint> a);
        void serializeBF(std::string BF_file);
        void deserializeBF(std::string BF_file);

        // void serialize1(std::string BF_file);

        int n;
        float p;
        int R;
        int k;
        // std::vector<bool> m_bits;
        // std::bitset<capacity> m_bits;
        bitArray* m_bits;
};

#endif
