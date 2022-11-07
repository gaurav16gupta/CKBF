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
#include <cstdlib>

#include <omp.h>
#include <filesystem>
#include <climits>
using namespace std;
#define  NUM_THREADS 1 // for yogi. run on 0-21,44-65

void print_array(vector <uint> arr) {
    for (int i=0;i < arr.size() ; i++) {
      cout << arr[i] << " ";
    }
    cout << endl;
}

void print_array(uint *  arr, int size) {
    for (int i=0;i <size ; i++) {
      cout << arr[i] << " ";
    }
    cout << endl;
}
//TODO(aditya) move to a config file
class Config {
        public:
                int range = 150000;
                int k = 5;
                string data_path = "/home/apd10/CKBF/data/";
                string results_path = "/home/apd10/CKBF/results/";
                //string file_id = "SRR649944";
                //string file_id = "SRR649955";
                string file_id = "SAMPLE";
                bool index = true;
                bool query = true;
                string bf_file;
                string data_file;
                string true_test_file;
                vector<string> false_test_files = {"random.query", "SAMPLE.query.p1", "SAMPLE.query.p2", "SAMPLE.query.p4"};

                // derived
                Config() {
                        bf_file = results_path + "BF_" + file_id + "_range_" + to_string(range) + "_k_" + to_string(k) + ".bf";
                        data_file = file_id + ".fastq";
                        true_test_file = "test_" + file_id + ".true.txt";
                        //false_test_file = "test_" + file_id + ".false.txt";
                }
                void print_config() {
                    cout << "RANGE: " << range << endl;
                    cout << "K: " << k << endl;
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

        //vector<std::string> keys = getFastqdata(mainfile);
        vector<string> keys = readlines(mainfile, 0); // The new format
        cout << "Total number of keys" << keys.size() << endl;
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



void query(BloomFilter * mybf, Config config, string test_file) {
        cout << "testing file " << test_file << endl;
        mybf->deserializeBF(config.bf_file);
        std::vector<string> testKeys = readlines(test_file, 0);
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
                                //break; // removing this so that we have fair latency comparison
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

void insert_ls_brute(BloomFilter * mybf, Config config, int token_size, int num_tokens, int local_size) {
        string mainfile = config.data_path + config.data_file;
        omp_set_num_threads(NUM_THREADS);
        chrono::time_point<chrono::high_resolution_clock> insert_t1 = chrono::high_resolution_clock::now();
        //vector<std::string> keys = getFastqdata(mainfile);
        vector<string> keys = readlines(mainfile, 0);
        cout << "Total number of keys" << keys.size() << endl;
        if (keys.size()==0){
                cout << mainfile << " does not exists or empty "<<endl;
                return;
        }
        int big_token_size = token_size + num_tokens - 1;
        #pragma omp parallel
        {
                int id = omp_get_thread_num();
                int nthd = omp_get_num_threads();
                int total_size = keys.size();
                int chunk = (total_size + nthd - 1) / nthd; 
                uint random_hash[config.k];
                uint min_hash[config.k];
                uint temp_hash[config.k];
                fill(min_hash, min_hash + config.k, INT_MAX);
                vector<uint> temp(config.k);

                print_array(min_hash, config.k);
                for( int i=id * chunk; (i < total_size) && (i < (id+1) * chunk); i++){
                        const char * key = keys[i].c_str();
                        //cout << "input " << key << endl;

                        fill(min_hash, min_hash + config.k, INT_MAX);
                        for (uint bi =0; bi < keys[i].size() - big_token_size + 1; bi++){
                                murmurhash(&key[bi], big_token_size, config.k, local_size, random_hash);
                                for (int li = 0; li < num_tokens; li++) {
                                    murmurhash(&key[bi + li], token_size, config.k, 1<<31, temp_hash);
                                    for (int j = 0; j < config.k; j++) {
                                        min_hash[j] = min(min_hash[j], temp_hash[j]);
                                    }
                                }
                                temp.clear();
                                for(int j=0; j<config.k;j++) {
                                    temp.push_back(min_hash[j] % (config.range - local_size + 1) + random_hash[j]);
                                }
                                //print_array(random_hash, config.k);
                                //print_array(min_hash, config.k);
                                //print_array(temp);
                                mybf->insert(temp);
                        }
                }
        }

        chrono::time_point<chrono::high_resolution_clock> insert_t2 = chrono::high_resolution_clock::now();
        cout << "Insertion time: " << time_seconds(insert_t1, insert_t2) << endl;
        mybf->serializeBF(config.bf_file);
        cout << "Packing in bloom filter: "<< mybf->m_bits->getcount()<< "/" << config.range << endl;
}

void query_ls_brute(BloomFilter * mybf, Config config, string test_file, int token_size, int num_tokens, int local_size) {
        cout << "testing file " << test_file << endl;
        mybf->deserializeBF(config.bf_file);
        std::vector<string> testKeys = readlines(test_file, 0);
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

        int big_token_size = token_size + num_tokens - 1;

        uint random_hash[config.k];
        uint min_hash[config.k];
        uint temp_hash[config.k];
        fill(min_hash, min_hash + config.k, INT_MAX);

        vector<uint> temp(config.k);
        for (std::size_t i=0; i<testKeys.size(); i++){
                membership = true;
                const char * key = testKeys[i].c_str();
                fill(min_hash, min_hash + config.k, INT_MAX);
                for (uint bi =0; bi < testKeys[i].size() - big_token_size + 1; bi++){
                        murmurhash(&key[bi], big_token_size, config.k, local_size, random_hash);
                        for (int li = 0; li < num_tokens; li++) {
                                murmurhash(&key[bi + li], token_size, config.k, 1<<31, temp_hash);
                                for (int j = 0; j < config.k; j++) {
                                        min_hash[j] = min(min_hash[j], temp_hash[j]);
                                }
                        }
                        // TODO(aditya) this should change. But using this as it is vector interface with bloomfilters
                        temp.clear();
                        for(int j=0; j<config.k;j++) {
                                temp.push_back(min_hash[j] % (config.range - local_size + 1) + random_hash[j]);
                        }

                        chrono::time_point<chrono::high_resolution_clock> check_t1 = chrono::high_resolution_clock::now();
                        this_membership = mybf->test(temp);
                        chrono::time_point<chrono::high_resolution_clock> check_t2 = chrono::high_resolution_clock::now();
                        check_time_accu += chrono::duration_cast<chrono::nanoseconds>(check_t2 - check_t1).count();
                        if (not this_membership){
                                membership = false;
                                //break;
                        }
                }
                if (membership){
                        pos = pos + 1;
                } 
                //cout << membership << " " << key << endl;
       }

        cout << "total check time is: " << check_time_accu << endl;
        cout << "observed positive rate: " << static_cast<float>(pos) / testKeys.size() << endl;
        chrono::time_point<chrono::high_resolution_clock> query_t2 = chrono::high_resolution_clock::now();
        cout <<"query time wall clock is :" << time_seconds(query_t1, query_t2) << endl;
}

int main(int argc, char** argv){
        int range, k;
        if (argc < 2) {
            cout << "Usage: " << "<program> range k [Need to move config to a file. but for temporary ease]"  << endl;
            return 0;
        }
        
        range = atoi(argv[1]);
        k = atoi(argv[2]);
        Config config;
        config.range = range;
        config.k = k;
        config.print_config();

        BloomFilter * mybf =  new BloomFilter(config.range, config.k);

        cout << endl << "===========random hash results =============" << endl;
        fflush(stdout);
        if (config.index){
                insert(mybf, config);
        }
        if (query){
                cout << "----- testing true strings ----" << endl;
                query(mybf, config, config.data_path + config.true_test_file);

                cout << "----- testing false strings ----" << endl;
                for (int i=0; i < config.false_test_files.size(); i++) {
                    query(mybf, config, config.data_path + config.false_test_files[i]);
                }
        }
        delete mybf;

        mybf =  new BloomFilter(config.range, config.k);

        cout << endl << "===========FuzzyLocality hash results =============" << endl;
        fflush(stdout);
        if (config.index){
                insert_ls_brute(mybf, config, 16, 16, 1024);
        }
        if (query){
                cout << "----- testing true strings ----" << endl;
                query_ls_brute(mybf, config, config.data_path + config.true_test_file, 16, 16, 1024);

                cout << "----- testing false strings ----" << endl;
                for (int i=0; i < config.false_test_files.size(); i++) {
                    query_ls_brute(mybf, config, config.data_path + config.false_test_files[i], 16, 16, 1024);
                }
        }

        return 0;
}
