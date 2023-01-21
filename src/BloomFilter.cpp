#include "BloomFilter.h"
#include <fcntl.h>
#include <sys/mman.h>

using namespace std;

BloomFilter::BloomFilter(uint64_t sz, uint32_t k_, bool disk, string name)
  : size(sz), k(k_) {
  if (disk) {
    file_write = open(("./results/"+name+".dat").c_str(), O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP); //todo: make sure results folder is preesent
    if (file_write<=0) cerr<<"can't open file "<<"./results/"+name+".dat"<<" to save BF on disk"<<endl;
    posix_fallocate(file_write, 0, sz >> 3); // sz >> 3 shd be multiple of 4096
    bits = reinterpret_cast<uint8_t*>(mmap(NULL, sz >> 3, PROT_WRITE, MAP_SHARED, file_write, 0));
  } else {
    bits = new uint8_t[sz >> 3];
  }
}


void BloomFilter::insert(uint64_t *hashes) {
  // TODO: hardcode for loop
  for (uint32_t i = 0; i < k; ++i) {
    bits[hashes[i] >> 3] |= 1 << (hashes[i] & 7);
  }
}

bool BloomFilter::test(uint64_t *hashes) {
  // TODO: compare with using if; hardcode for loop
  bool result = true;
  for (uint32_t i = 0; i < k; ++i) {
    result = result && bits[hashes[i] >> 3] & (1 << (hashes[i] & 7));
  }
  return result;
}

// void BloomFilter::release() {
//   munmap (bits, sz >> 3);
//   close (file_write);
// }

// make sure bits is all zeros
// https://bertvandenbroucke.netlify.app/2019/12/08/memory-mapping-files/

uint64_t BloomFilter::count() const {
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
