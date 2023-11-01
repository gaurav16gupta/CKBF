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
#include "MyBloom.h"
#include "MurmurHash3.h"
#include "Rambo_construction.h"
#include "utils.h"
#include "constants.h"
#include "bitArray.h"
#include <set>
#include <iterator>
#include <bitset>

using namespace std;

vector<uint> RAMBO:: hashfunc( std::string key, int len){
  // int hashvals[k];
  vector <uint> hashvals;
  uint op;
  for (int i=0; i<R; i++){
    MurmurHash3_x86_32(key.c_str(), len, 10, &op); //seed i
    hashvals.push_back(op%B);
  }
  return hashvals;
}


RAMBO::RAMBO(int n, int r1, int b1, int K){
  R = r1;
  B = b1;
  K = K;

  //range = ceil(-(n*log(p))/(log(2)*log(2))); //range
  range = n; //size of bloom filter (BFU)
  k = 2; //number of hash
  std::cout << "range" <<range<< '\n';

  Rambo_array = new BloomFiler*[B*R]; //array of pointers
  metaRambo = new vector<int>[B*R]; //constains set info in it.
  for(int b=0; b<B; b++){
    for(int r=0; r<R; r++){
      Rambo_array[b + B*r] = new BloomFiler(range, p, k);
    }
  }
}
// one time process- a preprocess step
void RAMBO::createMetaRambo(int K, bool verbose){
  for(int i=0;i<K;i++){
    vector<uint> hashvals = RAMBO::hashfunc(std::to_string(i), std::to_string(i).size()); // R hashvals, each with max value B
    for(int r=0; r<R; r++){
      metaRambo[hashvals[r] + B*r].push_back(i);
      // this->idx_and_r_to_b[i].push_back(hashvals[r]); // Assumes everything is in order
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
  #pragma omp parallel for num_threads(80)
  for(std::size_t i=0; i<keys.size(); ++i){
    for(int r=0; r<R; r++){
      vector<uint> temp = myhash(keys[i].c_str(), keys[i].size() , k,r, range);
      Rambo_array[hashvals[r] + B*r]->insert(temp);
    }
  }
}

// given inverted index type arrangement, kmer;files;files;..
void RAMBO::insertion2 (std::vector<string> alllines){
  //make this loop parallel
  //#pragma omp parallel for
  for(std::size_t i=0; i<alllines.size(); ++i){
    char d = ';';
    std::vector<string>KeySets =  line2array(alllines[i], d);//sets for a key

    std::vector<string>KeySet = line2array(KeySets[1], ',');
    for (uint j = 0; j<KeySet.size(); j++){
      vector<uint> hashvals = RAMBO::hashfunc(KeySet[j], KeySet[j].size()); // R hashvals
      for(int r=0; r<R; r++){
	      vector<uint> temp = myhash(KeySets[0].c_str(), KeySets[0].size() , k, r, range);// i is the key
        Rambo_array[hashvals[r] + B*r]->insert(temp);
      }
    }
  }
}

// give set and keys in the set
void RAMBO::insertionwithRead (std::string setID, std::string filenameSet){
  vector<uint> hashvals = RAMBO::hashfunc(setID, setID.size()); // R hashvals

  ifstream cntfile (filenameSet);
  // std::vector <std::string> allKeys;
  char keys[3000][32];
  int totKmerscnt = 0;
  
  while ( cntfile.good() ){
      string line1;
      while( getline ( cntfile, line1 ) ){

        for(int r=0; r<R; r++){
          vector<uint> temp = myhash(line1.substr(0, 31).c_str(), 31 , k,r, range);
          Rambo_array[hashvals[r] + B*r]->insert(temp);
        }
        totKmerscnt++;
        // if (totKmerscnt<3000){
        //   strcpy(keys[totKmerscnt], line1.substr(0, 31).c_str());
        //   }
        
        // totKmerscnt++;
        // if (totKmerscnt>=3000){
        //       #pragma omp parallel for
        //       for(std::size_t i=0; i<totKmerscnt; ++i){
        //         for(int r=0; r<R; r++){
        //           vector<uint> temp = myhash(keys[i], 31 , k,r, range);
        //           Rambo_array[hashvals[r] + B*r]->insert(temp);
        //         }
        //       }
        //     totKmerscnt =0;
        //   }
        
      }
    }
    std::cout<<"here3"<<endl;
}

bitArray RAMBO::query (string query_key, int len){
  // set<int> resUnion[R]; //constains union results in it.
  bitArray bitarray_K(Ki);
  std::set<unsigned int> to_check_next;
  // bitset<Ki> bitarray_K;
  // set<int> res;
  float count=0.0;
  vector<uint> check;
  for(int r=0; r<R; r++){
    if (r ==0){
      check = myhash(query_key, len , k, r,range); //hash values correspondign to the keys
      bitArray bitarray_K1(Ki);
      for(int b=0; b<B; b++){
          if (Rambo_array[b + B*r]->test(check)){
            for (auto idx: metaRambo[b + B*r]){
              // to_check_next.insert(this->idx_and_r_to_b[idx][0]);
              to_check_next.insert(idx_file_to_b[idx][0]);
              bitarray_K1.SetBit(idx);
          }
        }
      }
      bitarray_K = bitarray_K1;
    }
    else{
      std::set<unsigned int> to_check = to_check_next;
      to_check_next.clear();
      check = myhash(query_key, len , k, r,range); //hash values correspondign to the keys
      bitArray bitarray_K1(Ki);
      for(auto b: to_check){
          if (Rambo_array[b + B*r]->test(check)){
            for (auto idx: metaRambo[b + B*r]){
               if (r < R - 1) {
                  //  to_check_next.insert(this->idx_and_r_to_b[idx][r]);
                   to_check_next.insert(idx_file_to_b[idx][r]);
                  }
              bitarray_K1.SetBit(idx);
          }
        }
      }
      bitarray_K.ANDop(bitarray_K1.A);
    }
  }
  vector<uint>().swap(check);
  return bitarray_K;
}

void RAMBO::serializeRAMBO(string dir){
  for(int b=0; b<B; b++){
    for(int r=0; r<R; r++){
      string br = dir + to_string(b) + "_" + to_string(r) + ".txt";
      Rambo_array[b + B*r]->serializeBF(br);
    }
  }
}

void RAMBO::deserializeRAMBO(vector<string> dir){
  for(int b=0; b<B; b++){
    for(int r=0; r<R; r++){
      vector<string> br;
     	for (uint j=0; j<dir.size(); j++){
	  br.push_back(dir[j] + to_string(b) + "_" + to_string(r) + ".txt");
	}
      Rambo_array[b + B*r]->deserializeBF(br);

    }
  }
}
