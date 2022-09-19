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
#include "utils.h"
#include "bitArray.h"
#include <ctime>

using namespace std;

int main(int argc, char** argv){
string job(argv[1]);

int range = 2000000000; //size of BF
int k = 2;
std::cout << "range" <<range<< '\n';
// constructor
mybf = new BloomFiler(range, k);

float fp_ops;
float ins_time;
float query_time;
string SerOpFile ="results/SerFASTQ_" + to_string(range)+'_'+ to_string(k)+'.txt';

/////////////////INSERT//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
string mainfile = "/data/example.fastq";
chrono::time_point<chrono::high_resolution_clock> t3 = chrono::high_resolution_clock::now();
vector<std::string> keys = getFastqdata(mainfile);
if (keys.size()==0){
    std::cout<<mainfile<<" does not exists or empty "<<std::endl;
}
else{
  for(std::size_t i=0; i<keys.size(); ++i){
    for (uint x =0; x<keys[i].size()-31 +1; x++){
      vector<uint> temp = myhash(keys[i].substr(x, 31).c_str(), 31 , k,r, range);
      mybf.insert(temp);}
  }
}

chrono::time_point<chrono::high_resolution_clock> t4 = chrono::high_resolution_clock::now();
cout << chrono::duration_cast<chrono::nanoseconds>(t4-t3).count()/1000000000.0 << "sec\n";
ins_time = (t4-t3).count()/1000000000.0;

cout<<"Serializing RAMBO at: "<<SerOpFile<<endl;
mybf.serializeBF(SerOpFile);
cout<<"packing: "<<mybf.m_bits->getcount()<<endl;
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


///////////////Query/////////////////////////////////////////////////////////////
string SerOpFile ="results/SerFASTQ_" + to_string(range)+'_'+ to_string(k)+'.txt';
cout<<"deser starting"<<endl;
myRambo.deserializeRAMBO(SerOpFile);
std::cout << "desealized!" << '\n';

std::vector<string> testKeys = readlines("data/testfile.txt", 0);
cout<<"total number of queries : "<<testKeys.size()<<endl;

float pos=0;
// std::ofstream FPtestFile;
// FPtestFile.open("FPtestFileToy.txt");
std::clock_t t5_cpu = std::clock();
chrono::time_point<chrono::high_resolution_clock> t5 = chrono::high_resolution_clock::now();

for (std::size_t i=0; i<testKeys.size(); i++){
  check = myhash(testKeys[i], testKeys[i].size() , k, r,range); //hash values correspondign to the keys
  mybf.test(check){
    pos = pos + 1;
  }
}
cout<<"total pos is: "<<pos; // false positives/(all negatives)
cout<<"fp rate is: "<<(pos-posgt)/testKeys.size(); // false positives/(all negatives)

// cout<<endl;
std::clock_t t6_cpu = std::clock();
chrono::time_point<chrono::high_resolution_clock> t6 = chrono::high_resolution_clock::now();
float QTpt_cpu = 1000.0 * (t6_cpu-t5_cpu)/(CLOCKS_PER_SEC*testKeys.size()); //in ms
float QTpt = chrono::duration_cast<chrono::nanoseconds>(t6-t5).count()/(1000000.0*testKeys.size());
cout <<"query time wall clock is :" <<QTpt <<", cpu is :" <<QTpt_cpu<< " milisec\n\n";
// FPtestFile<<"query time wall clock is :" <<QTpt <<", cpu is :" <<QTpt_cpu<< " milisec\n\n";
// query_time = QTpt_cpu;
/////////////////////////////////////////////////////////////////////////////////
return 0;
}
