#include <iomanip>
#include <fstream>
#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>
#include <algorithm>
#include "IDLBloomFilter.h"
#include "Hasher.h"
#include "MurmurHash3.h"
#include "Rambo_constructionIDL.h"
#include "utils.h"
#include "constants.h"
#include "bitArray.h"
#include <set>
#include <iterator>
#include <bitset>
#include <omp.h>

using namespace std;

vector<uint> RAMBO:: hashfunc( std::string key, int len){
  // int hashvals[k];
  vector <uint> hashvals;
  uint op;
  for (int i=0; i<R; i++){
    MurmurHash3_x86_32(key.c_str(), len, i, &op); //seed i
    hashvals.push_back(op%B);
  }
  return hashvals;
}

RAMBO::RAMBO(int n, int r1, int b1, int K, int numThreads, string hashType, int universalHashRange){
  R = r1;
  B = b1;
  K = K;

  //range = ceil(-(n*log(p))/(log(2)*log(2))); //range
  range = n; //size of bloom filter (BFU)
  k = 4; //number of hash
  int seed =42;
  hasher = new Hasher*[numThreads];
  omp_set_num_threads(numThreads);
  if (hashType=="murmur"){
    for (uint32_t i = 0; i < numThreads; ++i) {
      hasher[i] = static_cast<Hasher*>(new MurmurHasher(range, k, seed));
    }
  }
  else if (hashType=="IDL"){
    for (uint32_t i = 0; i < numThreads; ++i) {
      int kMer=16;
      hasher[i] = static_cast<Hasher*>(new EfficientFuzzyHasherEXP(range, k, kMer, universalHashRange, seed));
    }
  }
  else{cout <<"wrong hash type!"<<endl;}

  Rambo_array = new IDLBloomFilter*[B*R]; //array of pointers
  metaRambo = new vector<int>[B*R]; //constains set info in it.
}

void RAMBO::initBFs(){
  for(int b=0; b<B; b++){
    for(int r=0; r<R; r++){
      Rambo_array[b + B*r] = new IDLBloomFilter(range, k, 0, 0, "temp");
    }
  }
}

// one time process- a preprocess step
void RAMBO::createMetaRambo(vector<string> setIDs, int K, bool verbose){
  for(int i=0;i<K;i++){
    vector<uint> hashvals = RAMBO::hashfunc(setIDs[i], setIDs[i].size()); // R hashvals, each with max value B
    for(int r=0; r<R; r++){
      metaRambo[hashvals[r] + B*r].push_back(i);
      idx_file_to_b[i][r] = hashvals[r];
    }
  }

  //print RAMBO meta deta
  if (verbose){
    for(int b=0; b<B; b++){
      for(int r=0; r<R; r++){
        for (auto it=metaRambo[b + B*r].begin(); it != metaRambo[b + B*r].end(); ++it)
          {std::cout << " " << *it;}
        {std::cout << "////";}
      }
      std::cout << '\n';
    }
}
}

// give set and keys in the set
void RAMBO::insertion (std::string setID, std::vector<std::string> keys){
  vector<uint> hashvals = RAMBO::hashfunc(setID, setID.size()); // R hashvals
  //make this loop parallel
  
  #pragma omp parallel for
  for(std::size_t i=0; i<keys.size(); ++i){
      // vector<uint> temp = myhash(keys[i].c_str(), keys[i].size() , k,r, range);
      uint64_t hashes[k];
      uint32_t threadId = omp_get_thread_num();
      // uint32_t threadId = 0;
      hasher[threadId]->setSequence(keys[i]);
      while (hasher[threadId]->hasNext()) {
          hasher[threadId]->hash(hashes);
          for(int r=0; r<R; r++){
            Rambo_array[hashvals[r] + B*r]->insert(hashes);
          }     
        }
    }
}

// // given inverted index type arrangement, kmer;files;files;..
// void RAMBO::insertion2 (std::vector<string> alllines){
//   //make this loop parallel
//   //#pragma omp parallel for
//   for(std::size_t i=0; i<alllines.size(); ++i){
//     char d = ';';
//     std::vector<string>KeySets =  line2array(alllines[i], d);//sets for a key

//     std::vector<string>KeySet = line2array(KeySets[1], ',');
//     for (uint j = 0; j<KeySet.size(); j++){
//       vector<uint> hashvals = RAMBO::hashfunc(KeySet[j], KeySet[j].size()); // R hashvals
//       for(int r=0; r<R; r++){
// 	      // vector<uint> temp = myhash(KeySets[0].c_str(), KeySets[0].size() , k, r, range);// i is the key
//         uint64_t hashes[k];
//         uint32_t threadId = 0;
//         hasher[threadId]->setSequence(keys[i]);
//         hasher[threadId]->hash(hashes);
//         Rambo_array[hashvals[r] + B*r]->insert(hashes);
//       }
//     }
//   }
// }

// // give set and keys in the set
// void RAMBO::insertionwithRead (std::string setID, std::string filenameSet){
//   vector<uint> hashvals = RAMBO::hashfunc(setID, setID.size()); // R hashvals

//   ifstream cntfile (filenameSet);
//   // std::vector <std::string> allKeys;
//   char keys[3000][32];
//   int totKmerscnt = 0;
  
//   while ( cntfile.good() ){
//       string line1;
//       while( getline ( cntfile, line1 ) ){

//         for(int r=0; r<R; r++){
//           // vector<uint> temp = myhash(line1.substr(0, 31).c_str(), 31 , k,r, range);
//           uint64_t hashes[k];
//           uint32_t threadId = 0;
//           hasher[threadId]->setSequence(keys[i]);
//           hasher[threadId]->hash(hashes);
//           Rambo_array[hashvals[r] + B*r]->insert(hashes);
//         }
//         totKmerscnt++;
//         // if (totKmerscnt<3000){
//         //   strcpy(keys[totKmerscnt], line1.substr(0, 31).c_str());
//         //   }
        
//         // totKmerscnt++;
//         // if (totKmerscnt>=3000){
//         //       #pragma omp parallel for
//         //       for(std::size_t i=0; i<totKmerscnt; ++i){
//         //         for(int r=0; r<R; r++){
//         //           vector<uint> temp = myhash(keys[i], 31 , k,r, range);
//         //           Rambo_array[hashvals[r] + B*r]->insert(temp);
//         //         }
//         //       }
//         //     totKmerscnt =0;
//         //   }
        
//       }
//     }
//     std::cout<<"here3"<<endl;
// }

bitArray RAMBO::query (string query_key, int len){
  // set<int> resUnion[R]; //constains union results in it.
  bitArray bitarray_K(Ki);
  std::set<unsigned int> to_check_next;
  // bitset<Ki> bitarray_K;
  // set<int> res;
  float count=0.0;
  vector<uint> check;
  uint32_t threadId = 0;
  hasher[threadId]->setSequence(query_key);
  uint32_t nHashes = hasher[threadId]->nHashesInSequence();
  uint64_t hashes[k * nHashes];

  for (uint32_t l = 0; l < nHashes; ++l) {
    hasher[threadId]->hash(hashes + l * k);
  }

//  int coutPos = 0;
  for(int r=0; r<R; r++){
      // check = myhash(query_key, len , k, r,range); //hash values correspondign to the keys
      bitArray bitarray_K1(Ki);
      for(int b=0; b<B; b++){
        bool positive = true;
        for (uint32_t l = 0; l < nHashes; ++l) {
            positive = positive && Rambo_array[b + B*r]->test(hashes + l * k, range);
        }
        if (positive){
          for (auto idx: metaRambo[b + B*r]){
            bitarray_K1.SetBit(idx);
          }
          // coutPos++;
        }
      }
      // cout<<coutPos<<" "; 
      if(r==0) bitarray_K = bitarray_K1;
      else bitarray_K.ANDop(bitarray_K1.A);
  }
  // cout<<endl; 
  vector<uint>().swap(check);
  return bitarray_K;
}
// if disk =1 save as .dat
void RAMBO::serializeRAMBO(string dir){
  for(int b=0; b<B; b++){
    for(int r=0; r<R; r++){
      string br = dir + to_string(b) + "_" + to_string(r) + ".dat";
      Rambo_array[b + B*r]->serializeBF(br);
    }
  }
}

void RAMBO::deserializeRAMBO(string dir){
  for(int b=0; b<B; b++){
    for(int r=0; r<R; r++){
      string br = dir + to_string(b) + "_" + to_string(r) + ".dat";
      Rambo_array[b + B*r] = new IDLBloomFilter(range, k, 0, 1, br);
      Rambo_array[b + B*r]->deserializeBF(br, 0);
    }
  }
}
