#include "MurmurHash3.h"
#include <iostream>
#include <cstring>
#include <chrono>
#include <vector>
#include "MyBloom.h"
#include <math.h>
#include <bitset>
#include "bitArray.h"
#include <algorithm> 

using namespace std;


BloomFiler::BloomFiler(int sz, int k){
      k = k; //number of hash
      m_bits = new bitArray(sz);
      }

void BloomFiler::insert(vector<uint> a){
  int N = a.size();
  for (int n =0 ; n<N; n++){
    m_bits->SetBit(a[n]);
  }
}

bool BloomFiler::test(vector<uint> a){
  int N = a.size();
  for (int n =0 ; n<N; n++){
      if (!m_bits->TestBit(a[n])){
        return false;
      }
  }
  return true;
}

void BloomFiler::serializeBF(string BF_file){
  m_bits->serializeBitAr(BF_file);
}

void BloomFiler::deserializeBF(string BF_file){
  m_bits->deserializeBitAr(BF_file);
}
