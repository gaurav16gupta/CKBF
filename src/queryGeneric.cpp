#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>
#include <algorithm>
#include <chrono>
#include "MyBloom.h"
#include "MurmurHash3.h"
#include "Rambo_construction.h"
#include "utils.h"
#include "constants.h"
#include "bitArray.h"
#include <ctime>

using namespace std;

int main(int argc, char** argv){

string job = "0";

int n_perSet = 2000000000; //cardinality of each set
int R_all = Ri;
int B_all = Bi;
int K = Ki; // total number of sets

// int R_all = atoi(argv[2]);
// int B_all = atoi(argv[3]);
// int K = atoi(argv[4]);

float fp_ops;
float ins_time;
float query_time;

// constructor
RAMBO myRambo(n_perSet, R_all, B_all, K);

//  details of RAMBO set partitioning
myRambo.createMetaRambo (K, false);

string SerOpFile ="results/RAMBO_Ser" + to_string(K)+'_'+job +'/';
vector<string> SerOpFile2;
SerOpFile2.push_back(SerOpFile); // mutliple files can be pushed here
myRambo.deserializeRAMBO(SerOpFile2);
std::cout << "desealized!" << '\n';

//get the ecoliNames 
std::vector<string> ecoliNames;
for (int batch =0; batch<K/100; batch++){
    string dataPath = "data/"+ job +"/" + to_string(batch) + "_indexed.txt";
    std::vector<string> alllines = readlines(dataPath, 0);
    for(uint i=0; i<alllines.size(); i++){
        ecoliNames.push_back(line2array(alllines[i], ',')[1]);  
    }
}
label1:
    // test RAMBO
    string queryfile;
    cout<<"Enter query or filename: ";
    cin >> queryfile;
    std::vector<string> testKeys;
    bool verbose=0;
    if(queryfile.substr(queryfile.find_last_of(".") + 1) == "txt"){
        testKeys = readlines(queryfile, 0);
        //std::vector<string> alllines = readlines(queryfile, 0);
        //for(uint i=0; i<alllines.size(); i++){
        //    testKeys.push_back(line2array(alllines[i], ' ')[0]);
        //}
        cout<<"total number of queries : "<<testKeys.size()<<endl;
    }
    else{
        testKeys.push_back(queryfile);
        verbose =1;
    }

    // float fp=0;
    std::ofstream FPtestFile; //log file
    FPtestFile.open("FPtestFileToy.txt");
    std::clock_t t5_cpu = std::clock();
    chrono::time_point<chrono::high_resolution_clock> t5 = chrono::high_resolution_clock::now();
    vector <uint> ops;
    cout<<testKeys.size()<<endl;
    for (std::size_t i=0; i<testKeys.size(); i++){
            //cout<<testKeys[i]<<endl;
        bitArray MemVec = myRambo.query(testKeys[i], testKeys[i].size());
        ops = MemVec.get1locs();
        if (verbose){
            for (std::size_t j=0; j<ops.size(); j++){
                cout<<ecoliNames[ops[j]]<<".out ";
            }}
        // cout<<ops.size()<<endl;
        // break;
        // if (ops.size()<7){cout<<testKeys[i]<<endl;}
        }
    cout<<endl;

    std::clock_t t6_cpu = std::clock();
    chrono::time_point<chrono::high_resolution_clock> t6 = chrono::high_resolution_clock::now();

    float QTpt_cpu = 1000.0 * (t6_cpu-t5_cpu)/(CLOCKS_PER_SEC); //in ms
    float QTpt = chrono::duration_cast<chrono::nanoseconds>(t6-t5).count()/(1000000.0);
    cout <<"query time is :" <<QTpt/1000.0 <<" sec for "<<to_string(testKeys.size())+ " queries = "<< QTpt/testKeys.size()<<" millisec per query \n\n";
    FPtestFile<<"query time is :" <<QTpt <<" milisec\n\n";
    query_time = QTpt_cpu;
goto label1;

return 0;
}

