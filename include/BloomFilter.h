#ifndef _MYBLOOM_
#define _MYBLOOM_

#include <iostream>

class BloomFilter{
public:
    BloomFilter(uint32_t sz, uint32_t k);
    void insert(uint32_t *hashes);
    bool test(uint32_t *hashes);

    uint8_t *bits;
    uint32_t k;
};

#endif
