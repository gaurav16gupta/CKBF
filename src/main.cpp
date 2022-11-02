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
#include <omp.h>
#include <filesystem>

using namespace std;
#define  NUM_THREADS 44 // for yogi. run on 0-21,44-65

//TODO(aditya) move to a config file
class Config {
    public:
        int range = 2000000000;
        int k = 7;
        string data_path = "/home/apd10/CKBF/data/";
        string results_path = "/home/apd10/CKBF/results/";
        string file_id = "SRR649944";
        bool index = true;
        bool query = true;
        string bf_file;
        string data_file;
        string test_file;

        // derived
        Config() {      
              bf_file = results_path + "BF_" + file_id + "_range_" + to_string(range) + "_k_" + to_string(k) + ".bf";
              data_file = file_id + ".fastq";
              test_file = "test_" + file_id + ".txt";
        }
};

float time_seconds(chrono::time_point<chrono::high_resolution_clock> start,
                   chrono::time_point<chrono::high_resolution_clock> end) {

        float time = chrono::duration_cast<chrono::nanoseconds>(end - start).count()/1000000000.0 ;
        return time;
}

void insert(BloomFilter * mybf, Config config) {
        string mainfile = config.data_path + config.data_file;
        omp_set_num_threads(NUM_THREADS);
        chrono::time_point<chrono::high_resolution_clock> insert_t1 = chrono::high_resolution_clock::now();

        vector<std::string> keys = getFastqdata(mainfile);
        cout << "Total number of keys" << keys.size() << endl;
        cout << keys[0] <<endl;
        if (keys.size()==0){
                cout << mainfile << " does not exists or empty "<<endl;
                return;
        }
        #pragma omp parallel
        {
                int i;
                int id = omp_get_thread_num();
                int nthd = omp_get_num_threads();
                int total_size = keys.size();
                int chunk = (total_size + nthd - 1) / nthd; 
                // cout<<nthd<<endl;
                for( i=id * chunk; (i < total_size) && (i < (id+1) * chunk); i++){
                        for (uint x =0; x<keys[i].size()-31 +1; x++){
                                // vector<uint> temp = xxhash32(keys[i].substr(x, 31).c_str(), 31 , k, range);
                                vector<uint> temp = myhash2(keys[i].substr(x, 31).c_str(), 31, config.k, config.range);
                                mybf->insert(temp);
                        }
                }
        }

        chrono::time_point<chrono::high_resolution_clock> insert_t2 = chrono::high_resolution_clock::now();
        cout << "Insertion time: " << time_seconds(insert_t1, insert_t2) << endl;
        mybf->serializeBF(config.bf_file);
        cout << "Packing in bloom filter: "<< mybf->m_bits->getcount()<< "/" << config.range << endl;
}

void query(BloomFilter * mybf, Config config) {
                mybf->deserializeBF(config.bf_file);
                std::vector<string> testKeys = readlines(config.data_path + config.test_file, 0);
                cout << "total number of queries : "<<testKeys.size()<<endl;
                cout << testKeys[0] <<endl;

                float pos=0;
                chrono::time_point<chrono::high_resolution_clock> query_t1 = chrono::high_resolution_clock::now();

                vector<uint> check;
                bool membership;
                bool this_membership;
                float hash=0;
                float look=0;
                int64_t hash_time_accu = 0;
                int64_t check_time_accu = 0;

                for (std::size_t i=0; i<testKeys.size(); i++){
                        membership = true;
                        for (uint q =0; q < testKeys[i].size()-31 +1; q++){

                                chrono::time_point<chrono::high_resolution_clock> hash_t1 = chrono::high_resolution_clock::now();
                                check = myhash2(testKeys[i].substr(q, 31), 31, config.k, config.range);
                                chrono::time_point<chrono::high_resolution_clock> hash_t2 = chrono::high_resolution_clock::now();
                                hash_time_accu += chrono::duration_cast<chrono::nanoseconds>(hash_t2 - hash_t1).count();

                                chrono::time_point<chrono::high_resolution_clock> check_t1 = chrono::high_resolution_clock::now();
                                this_membership = mybf->test(check);
                                chrono::time_point<chrono::high_resolution_clock> check_t2 = chrono::high_resolution_clock::now();
                                check_time_accu += chrono::duration_cast<chrono::nanoseconds>(check_t2 - check_t1).count();
                                if (not this_membership){
                                        membership = false;
                                        break;
                                }
                        }
                        if (membership){
                                pos = pos + 1;
                        }
                }
                cout << "total hash time is: " << hash_time_accu << endl;
                cout << "total check time is: " << check_time_accu << endl;
                cout << "observed positive rate: " << static_cast<float>(pos) / testKeys.size() << endl;
                chrono::time_point<chrono::high_resolution_clock> query_t2 = chrono::high_resolution_clock::now();
                cout <<"query time wall clock is :" << time_seconds(query_t1, query_t2) << endl;
}

int main(int argc, char** argv){
        Config config;
        BloomFilter * mybf =  new BloomFilter(config.range, config.k);
        
        if (config.index){
            insert(mybf, config);
        }
        if (query){
            query(mybf, config);
        }
        return 0;
}

