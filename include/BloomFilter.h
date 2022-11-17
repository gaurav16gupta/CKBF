#ifndef _MYBLOOM_
#define _MYBLOOM_
#include <atomic>
#include <iostream>
#include <fcntl.h>

class BloomFilter{
public:
    BloomFilter(uint64_t sz, uint32_t k_, bool disk);
    void insert(uint32_t *hashes);
    bool test(uint32_t *hashes);
    uint32_t count() const;
    // void release();

    uint8_t *bits;
    int file_write;
    uint32_t size;
    uint32_t k;
};

#endif



