#include "BloomFilter.h"

using namespace std;

BloomFilter::BloomFilter(uint32_t sz, uint32_t k)
  : bits(new uint8_t[sz >> 3]), k(k) {

}

void BloomFilter::insert(uint32_t *hashes) {
  // TODO: hardcode for loop
  for (uint32_t i = 0; i < k; ++i) {
    bits[hashes[i] >> 3] |= 1 << (hashes[i] & 7);
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
