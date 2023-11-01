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

int n_perSet =  atoi(argv[3]); //Bloom filter (BFU) size, range
int R_all = Ri;
int B_all = Bi;

int K = Ki; // total number of sets

float fp_ops;
float ins_time;
float query_time;

// constructor
// RAMBO myRambo(n_perSet, R_all, B_all, K);
RAMBO myRambo(n_perSet, R_all, B_all, K, 1, algo, UhRange);

//  details of RAMBO set partitioning
string dataPath = "../CKBF/data/fastqFiles/names.txt";
std::vector<string> setIDs = readlines(dataPath, 0);
myRambo.createMetaRambo (setIDs, K, false);

// string SerOpFile ="results1BR2/RAMBO_Ser_k4_" +to_string(UhRange)+"_"+ to_string(K)+"_"+algo+'/';
string SerOpFile ="results/RAMBO_Ser_" + to_string(K)+"_"+algo+'_'+ to_string(n_perSet)+"_"+to_string(UhRange)+'/';
cout<<SerOpFile<<endl;


// (deser)
  cout<<"deser starting"<<endl;
  myRambo.deserializeRAMBO(SerOpFile);
  std::cout << "desealized!" << '\n';

// (test)
    // test RAMBO
    // std::vector<string> alllines = readlines("data/ArtfcKmersToy"+to_string(K)+".txt", 0);
    // std::vector<string> alllines = readlines("../CKBF/data/queries/DRR021939query.p1", 0);
    std::vector<string> alllines = readlines("data/queryIDL.txt", 0);

    std::vector<string> testKeys;
    std::vector<int> gt_size;
    for(uint i=0; i<alllines.size(); i++){
            // std::vector<string>KeySets =  line2array(alllines[i], ';');//sets for a key
            // testKeys.push_back(KeySets[0]);
            // gt_size.push_back( line2array(KeySets[1], ',').size() );
            testKeys.push_back(alllines[i]);
            gt_size.push_back( 0 );
    }
    // myRambo.createMetaRambo (K, false);
    // cout<<"load: "<<myRambo.Rambo_array[0]->m_bits->getcount();
    cout<<"total number of queries : "<<testKeys.size()<<endl;
        // myRambo.insertion2 (alllines);

	float fp=0;
	std::ofstream FPtestFile;
	FPtestFile.open("FPtestFileToy.txt");
	std::clock_t t5_cpu = std::clock();
	chrono::time_point<chrono::high_resolution_clock> t5 = chrono::high_resolution_clock::now();

	for (std::size_t i=0; i<testKeys.size(); i++){
			bitArray MemVec = myRambo.query(testKeys[i], testKeys[i].size());
			// cout<<MemVec.getcount()<<endl;
			// cout<<gt_size[i]<<endl;
			fp = fp + (MemVec.getcount())*0.1/((K - gt_size[i])*0.1);
		}

    std::clock_t t6_cpu = std::clock();
	chrono::time_point<chrono::high_resolution_clock> t6 = chrono::high_resolution_clock::now();

	cout<<"fp rate is: "<<fp/(testKeys.size()); // false positives/(all negatives)
	FPtestFile<<"fp rate is: "<<fp/(testKeys.size()); // false positives/(all negatives)

	fp_ops = fp/(testKeys.size());

	cout<<endl;
	
	float QTpt_cpu = 1000.0 * (t6_cpu-t5_cpu)/(CLOCKS_PER_SEC*testKeys.size()); //in ms
	float QTpt = chrono::duration_cast<chrono::nanoseconds>(t6-t5).count()/(1000000.0*testKeys.size());
	cout <<"query time wall clock (milisec)\n";
    cout <<QTpt <<endl;
  	cout<<"query time wall clock is :" <<QTpt <<", cpu is :" <<QTpt_cpu<< " milisec\n\n";

	FPtestFile<<"query time wall clock is :" <<QTpt <<", cpu is :" <<QTpt_cpu<< " milisec\n\n";
	query_time = QTpt_cpu;

return 0;
}
