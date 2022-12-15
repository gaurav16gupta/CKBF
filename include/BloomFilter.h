#ifndef _MYBLOOM_
#define _MYBLOOM_

#include <iostream>

class BloomFilter{
public:
    BloomFilter(uint64_t sz, uint32_t k);
    void insert(uint32_t *hashes);
    bool test(uint32_t *hashes);
    uint64_t count() const;

    uint8_t *bits;
    uint64_t size;
    uint32_t k;
};

#endif
