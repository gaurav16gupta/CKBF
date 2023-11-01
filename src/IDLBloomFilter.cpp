#include "IDLBloomFilter.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <sstream>

using namespace std;

IDLBloomFilter::IDLBloomFilter(uint64_t sz, uint32_t k_, bool disk, bool readOnly, string filePath)
  : size(sz), k(k_) {

  if (disk) {
    if (readOnly) {
      file_ = open(filePath.c_str(), O_RDONLY, 0);
      if (file_<=0) cerr<<"can't open file "<<filePath<<" to load BF on disk"<<endl;
      bits = reinterpret_cast<uint8_t*>(mmap(NULL, sz >> 3, PROT_READ, MAP_SHARED, file_, 0));
    } else {
      file_ = open(filePath.c_str(), O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP); //todo: make sure results folder is preesent
      if (file_<=0) cerr<<"can't open file "<<filePath<<" to save BF on disk"<<endl;
      posix_fallocate(file_, 0, sz >> 3); // sz >> 3 shd be multiple of 4096
      bits = reinterpret_cast<uint8_t*>(mmap(NULL, sz >> 3, PROT_WRITE, MAP_SHARED, file_, 0));
    }
  } else {
    bits = new uint8_t[sz >> 3];
    
  }

  if (!readOnly) {
    // initialize the bits to 0
    for (uint64_t i = 0; i < (sz>>3); ++i) {
      bits[i] = 0;
    }
  }
}


void IDLBloomFilter::insert(uint64_t *hashes) {
  // TODO: hardcode for loop
  for (uint32_t i = 0; i < k; ++i) {
    bits[hashes[i] >> 3] |= 1 << (hashes[i] & 7);
  }
}

bool IDLBloomFilter::test(uint64_t *hashes) {
  // TODO: compare with using if; hardcode for loop
  bool result = true;
  for (uint32_t i = 0; i < k; ++i) {
    result = result && bits[hashes[i] >> 3] & (1 << (hashes[i] & 7));
  }
  return result;
}

bool IDLBloomFilter::test(uint64_t *hashes, uint64_t mod) {
  // TODO: compare with using if; hardcode for loop
  bool result = true;
  
  for (uint32_t i = 0; i < k; ++i) {
    result = result && bits[(hashes[i] % mod) >> 3] & (1 << ((hashes[i] % mod) & 7));
  }
  return result;
}

void IDLBloomFilter::release() {
  munmap (bits, size >> 3);
}
// make sure bits is all zeros
// https://bertvandenbroucke.netlify.app/2019/12/08/memory-mapping-files/

uint64_t IDLBloomFilter::count() const {
  uint64_t cnt = 0;
  for (uint64_t i = 0; i < size >> 3; ++i) {
    uint8_t num = bits[i];
    while (num > 0) {
      cnt += num % 2;
      num /= 2;
    }
  }
  return cnt;
}

void IDLBloomFilter::serializeBF(string BF_file){
  file_ = open(BF_file.c_str(), O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP); //todo: make sure results folder is preesent
  if (file_<=0) cerr<<"can't open file "<<BF_file<<" to save BF on disk"<<endl;
  posix_fallocate(file_, 0, size >> 3);
  uint8_t *bitsDisk;
  bitsDisk = reinterpret_cast<uint8_t*>(mmap(NULL, size >> 3, PROT_WRITE, MAP_SHARED, file_, 0));
  //memcopy here instead
  for (uint64_t i = 0; i < size >> 3; ++i) {
    bitsDisk[i] =bits[i];
  }
}

void IDLBloomFilter::deserializeBF(string BF_file, bool disk){
  file_ = open(BF_file.c_str(), O_RDONLY, 0);
  if (file_<=0) cerr<<"can't open file "<<BF_file<<" to load BF on disk"<<endl;
  if (disk) 
    bits = reinterpret_cast<uint8_t*>(mmap(NULL, size >> 3, PROT_READ, MAP_SHARED, file_, 0));
  else{
    uint8_t *bitsDisk;
    bitsDisk = reinterpret_cast<uint8_t*>(mmap(NULL, size >> 3, PROT_READ, MAP_SHARED, file_, 0));
    //memcopy here instead
    for (uint64_t i = 0; i < size >> 3; ++i) {
      bits[i] =bitsDisk[i];
    }
  }
}
