#include "BloomFilter.h"

using namespace std;

BloomFilter::BloomFilter(uint32_t sz, uint32_t k)
  : bits(new uint8_t[sz >> 3]), size(sz), k(k) {

}

void BloomFilter::insert(uint32_t *hashes) {
  // TODO: hardcode for loop
  for (uint32_t i = 0; i < k; ++i) {
    bits[hashes[i] >> 3] |= 1 << (hashes[i] & 7);
    #pragma omp flush(bits)

    while ((bits[hashes[i] >> 3] & (1 << (hashes[i] & 7))) == 0) {
      cerr << "set bit again" << endl;
      bits[hashes[i] >> 3] |= 1 << (hashes[i] & 7);
      #pragma omp flush(bits)
    }
  }
}

bool BloomFilter::test(uint32_t *hashes) {
  // TODO: compare with using if; hardcode for loop
  bool result = true;
  for (uint32_t i = 0; i < k; ++i) {
    result = result && bits[hashes[i] >> 3] & (1 << (hashes[i] & 7));
  }
  return result;
}

uint32_t BloomFilter::count() const {
  uint32_t cnt = 0;
  for (uint32_t i = 0; i < size >> 3; ++i) {
    uint8_t num = bits[i];
    while (num > 0) {
      cnt += num % 2;
      num /= 2;
    }
  }
  return cnt;
}
