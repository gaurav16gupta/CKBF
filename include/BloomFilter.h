#ifndef _MYBLOOM_
#define _MYBLOOM_
#include <atomic>
#include <iostream>
#include <fcntl.h>

class BloomFilter{
public:
    BloomFilter(uint64_t sz, uint32_t k_, bool disk, bool readOnly, std::string filePath);
    void insert(uint64_t *hashes);
    bool test(uint64_t *hashes);
    bool test(uint64_t *hashes, uint64_t mod);
    void release();
    uint64_t count() const;
    // void release();

    uint8_t *bits;
    int file_;
    uint64_t size;
    uint32_t k;
};

#endif



