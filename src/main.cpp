#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>
#include <algorithm>
#include <chrono>
#include "IDLBloomFilter.h"
#include "Hasher.h"
#include "MurmurHash3.h"
#include "Rambo_constructionIDL.h"
#include "utils.h"
#include "constants.h"
#include "bitArray.h"
#include <ctime>

using namespace std;

int main(int argc, char** argv){

string algo(argv[1]);
int UhRange = atoi(argv[2]);

int n_perSet = atoi(argv[3]); // Bloom filter (BFU) size, range
int R_all = Ri;
int B_all = Bi;

int K = Ki; // total number of sets

float fp_ops;
float ins_time;
float query_time;

// constructor
// RAMBO myRambo(n_perSet, R_all, B_all, K);
RAMBO myRambo(n_perSet, R_all, B_all, K, 64, algo, UhRange);

//  details of RAMBO set partitioning
string dataPath = "../CKBF/data/fastqFiles/names.txt";
std::vector<string> setIDs = readlines(dataPath, 0);
myRambo.createMetaRambo (setIDs, K, false);

//insert itno RAMBO

// string SerOpFile ="results/RAMBO_Ser_k4_" +to_string(UhRange)+"_"+ to_string(K)+"_"+algo+'/';
string SerOpFile ="results/RAMBO_Ser_" + to_string(K)+"_"+algo+'_'+ to_string(n_perSet)+"_"+to_string(UhRange)+'/';

cout<<SerOpFile<<endl;

myRambo.initBFs();
// insert

  //log files
  std::ofstream failedFiles;
  failedFiles.open("logFileToy_"+ to_string(K)+".txt");
  int stpCnt = 0;

    chrono::time_point<chrono::high_resolution_clock> t3 = chrono::high_resolution_clock::now();    
    // #pragma omp parallel for num_threads(80)
    for (uint ss=0; ss<K; ss++){ //for every file
      stpCnt++;
      char d = ',';
      // vector<std::string> setID = line2array(setIDs[ss], d);
      string mainfile = "../CKBF/data/fastqFiles/" + setIDs[ss];

      // chrono::time_point<chrono::high_resolution_clock> cp1 = chrono::high_resolution_clock::now();
      // vector<std::string> keys = getctxdata(mainfile);
      vector<std::string> keys = getFastqdata(mainfile);
      
      myRambo.insertion(setIDs[ss], keys);
      // chrono::time_point<chrono::high_resolution_clock> cp4 = chrono::high_resolution_clock::now();

      failedFiles<<mainfile<<" : "<<keys.size()<<std::endl;
      if (keys.size()==0){
          std::cout<<mainfile<<" does not exists or empty "<<std::endl;
          failedFiles<<mainfile<<" does not exists or empty "<<std::endl;
      }
      // cout<<setIDs[ss]<<": "<<keys.size()<<endl;
      // cout <<"ctx read: "<< chrono::duration_cast<chrono::nanoseconds>(cp2-cp1).count()/1000000000.0 << "sec\n";
      // cout <<"insert: "<< chrono::duration_cast<chrono::nanoseconds>(cp4-cp1).count()/1000000000.0 << "sec\n";

    }
    chrono::time_point<chrono::high_resolution_clock> t4 = chrono::high_resolution_clock::now();
    cout << chrono::duration_cast<chrono::nanoseconds>(t4-t3).count()/1000000000.0 << "sec\n";
    ins_time = (t4-t3).count()/1000000000.0;
    failedFiles<<"insertion time (sec) of 100 files: "<<ins_time<<endl;

  //serialize myRambo
  try{
    string command = "mkdir "+ SerOpFile;
    system(command.c_str());
  }
  catch(int myNum){
    cout <<"folder present"<<endl;
  }

	cout<<"Serializing RAMBO at: "<<SerOpFile<<endl;
    myRambo.serializeRAMBO(SerOpFile);
   //gives number of 1s in 30 BFs
  //  for (int kpp=0;kpp<B_all*R_all;kpp++){
	// cout<<"packing: "<<kpp<<" : "<<myRambo.Rambo_array[kpp]->count()<<endl;
	//  }
  

return 0;
}
