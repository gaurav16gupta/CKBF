//contains all the other non-RAMBO functions
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>
#include <algorithm>
#include <chrono>
#include <random>
#include "MyBloom.h"
#include "MurmurHash3.h"
#include "xxhash32.h"
#include "utils.h"
#include <map>
#include <cassert>
using namespace std;


std::vector <std::string>  getFastqdata(string filenameSet){
  //get the size of Bloom filter by count
  ifstream cntfile (filenameSet);
  std::vector <std::string> allLines;
  int totKmerscnt = 0;
  while (cntfile.good()){
    string line1;
    int a =0;
    while( getline ( cntfile, line1 ) ){
    int sz = line1.size(); 
           
    if (a%4==1 && sz>=31){
      allLines.push_back(line1);
      totKmerscnt=totKmerscnt+sz-31;
      }
    a++;	   
    }
  }
  std::cout<<"total kmers in this file: "<<totKmerscnt<<std::endl;
  return allLines;
}


void murmurhash( const char *  key, int len, int k, uint range, uint * hashvals){
  uint op; // takes 4 byte
  for (int i=0; i<k; i++){
    MurmurHash3_x86_32(key, len, i + 10244231, &op);
    hashvals[i] = op%range;
  }
}

vector<uint> myhash( std::string key, int len, int k, int range){
  vector <uint> hashvals;
  // should be a pointer array
  uint op; // takes 4 byte
  for (int i=0; i<k; i++){
    MurmurHash3_x86_32(key.c_str(), len, i, &op);
    hashvals.push_back(op%range);
    // push_back might be slow due to memory reallocation
  }
  return hashvals;
}

vector<uint> myhash2( std::string &&key, int len, int k, int range){
  vector <uint> hashvals(k);
  uint op; // takes 4 byte
  for (int i=0; i<k; i++){
    MurmurHash3_x86_32(key.c_str(), len, i, &op);
    hashvals[i] = op%range;
  }
  return hashvals;
}

vector<uint> xxhash32(std::string &&key, int len, int k, int range) {
  vector<uint> hashvals(k);
  for (int i = 0; i < k; ++i) {
    hashvals[i] = XXHash32::hash(key.c_str(), len, i) % range;
  }
  return hashvals;
}

uint *rand_nums(int len, int range) {
  uint *rn = new uint[len];
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist(0, range-1); // distribution in range [1, 6]
  for (int i = 0; i < len; ++i) {
    rn[i] = dist(rng);
  }
  return rn;
}

vector<uint> tabulation_hash(std::string &&key, int k, uint rand_nums[]) {
  vector<uint> hashvals(k);
  for (int i = 0; i < k; ++i) {
    uint hv = 0;
    for (int j = 0; j < 31; ++j) {
      hv ^= rand_nums[256 * j + key[j]];
    }
    cout << hv << ' ';
    hashvals[i] = hv;
  }
  return hashvals;
}

vector<uint> minhash( std::string key, int len, int k, int range){
  vector <uint> hashvals;
  vector <uint> minhashvals;
  uint op; // takes 4 byte
  for (int i=0; i<k; i++){
    for (uint x =0; x<key.size()-10 +1; x++){
      MurmurHash3_x86_32(key.substr(x, 10).c_str(), 10, i, &op);
      minhashvals.push_back(op);
    }
    op = *min_element(minhashvals.begin(), minhashvals.end());
    hashvals.push_back(op%range);
  }
  return hashvals;
}

vector<uint> myhashCheap( std::string key, int len, int k, int range){
  vector <uint> hashvals;
  uint op; // takes 4 byte
  uint a[3] = {911,523,571};
  uint b[3] = {861719,332803,9397};
  for (int i=0; i<k; i++){
    op=0;
    for (int j=0; j<len; j++){
      op= op+ pow(3,j)*(a[i]*key[j]+ b[i]);
    }
    
    hashvals.push_back(op%range);
  }
  return hashvals;
}

//num control how many lines to get
std::vector<string> readlines( string path, int num){
  ifstream pathfile (path);
  std::vector <std::string> allfiles;
  int count=0;
  while ( pathfile.good() )
    {
      string line1;
      while( getline (pathfile, line1 ) ){
        count++;
        allfiles.push_back(line1);
        if (count >num && num){
          break;
        }
        }
  std::cout << count<< '\n';
  return allfiles;
  }
}


// //readlines from a file
// std::vector<string> getsets( string path){
//   //get the size of Bloom filter by count
//   ifstream cntfile (path);
//   std::vector <std::string> allKeys;
//   int linecnt = 0;
//   while (cntfile.good())
//       {
//           string line1, vals;
//           while( getline ( cntfile, line1 ) ){
//               stringstream is;
//               is<<line1;
//               if (linecnt==0 ){
//                 while(getline (is, vals, ' ' )){
//                   allKeys.push_back(vals);}
//               }
//               else{
//                 while(getline (is, vals, ' ' )){
//                   allKeys.push_back(vals);
//                   break;}
//               }
//                 linecnt++;
//           }
//       }
//       std::cout<<"total lines from this file: "<<linecnt-3<<std::endl;
//       return allKeys;
// }

// // parse string into array given delimiter
// std::vector<string> line2array( string line, char d){
//   stringstream is;
//   is<<line;
//   std::vector <std::string> op;
//   string vals;
//   while( getline (is, vals, d)){
//     op.push_back(vals);
//   }
//   return op;
// }


// //file write
// void writeRAMBOresults(string path, int rows, int cols, float* values){
//   ofstream myfile;
//   myfile.open (path);
//   for (int i =0;i<rows; i++){
//     for (int j =0;j<cols; j++){
//       myfile << to_string(values[i*cols +j])<<",";
//     }
//     myfile <<"\n";
//   }
//   myfile.close();
// }


// std::vector<string> getRandomTestKeys(int keysize, int n){
//   static const char alphanum[] = "ATGC";
//   std::vector<string> s;

//   for (int j = 0; j < n; ++j){
//     string st;
//       for (int i = 0; i < keysize; ++i) {
//           st = st + alphanum[rand()%4];
//       }
//       s.push_back(st);
//   }
// return s;
// }

// std::map<std::string, vector<int>> makeInvIndex(int n, vector<std::string> foldernames){
//   std::map<std::string, vector<int>> m;
//   for (uint f=0; f<foldernames.size(); f++){
//     string foldername = foldernames[f];
//     for (int batch=0; batch<47; batch++){
//   	string dataPath = foldername + to_string(batch)+ "_indexed.txt";
//   	std::vector<string> setIDs = readlines(dataPath, 0);
//   	cout<<setIDs[0]<<endl;

//   	for (uint ss=0; ss<setIDs.size(); ss++){
//     	  char d = ',';
//     	  vector<std::string> setID = line2array(setIDs[ss], d);
//           string mainfile = foldername + setID[1]+ ".out";
// 	  cout<<"getting keys"<<endl;
//     	  vector<std::string> keys = getctxdata(mainfile);
//     	  cout<<"gotkeys"<<endl;

// 	  if (ss==0 && batch ==0 && f ==0){
//       		for (int i =0; i<n; i++){
//         	   m[keys[i]].push_back(std::stoi(setID[0]));
//       		}
//       	        cout<<"completed first itr"<<endl;
//  		for(map<string, vector<int> >::iterator it = m.begin(); it != m.end(); ++it){
// 		   std::cout << it->first <<it->second[0]<<"\n";
// 		}
//     	   }
//     	  else{

//     		for (uint i =0; i<keys.size(); i++){
//       		   if (m.find(keys[i]) != m.end()){
//         	     std::cout << "map contains the key!\n";
//         	     m[keys[i]].push_back(std::stoi(setID[0]));
//       		   }
//     		}
// 	   }
//           cout<<ss<<endl;
//   	}
//      }
//   }
//   return m;
// }

// std::vector<std::string> getkmers(std::string query_key, int kmersize){
//   std::vector<std::string> query_kmers;
//   for (uint idx =0; idx<query_key.size()-kmersize +1; idx++){
//     query_kmers.push_back(query_key.substr(idx, kmersize));
//   }
//  return query_kmers;
// }

