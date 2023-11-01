#ifndef _RamboConstruction_
#define _RamboConstruction_
#include <vector>
#include <set>
#include <string>
#include <bitset>
#include "constants.h"
#include "bitArray.h"


// vector<uint> hashfunc( void* key, int len, int R, int B){
// }

class RAMBO{
    public:

        // RAMBO(int n, int r1, int b1, int K);
        RAMBO(int n, int r1, int b1, int K, int numThreads, string hashType, int universalHashRange);
        void initBFs();
        std::vector<uint> hashfunc( std::string key, int len);
        void insertion (std::string setID, std::vector<std::string> keys);
        std::set<int> takeunion(std::set<int> set1, std::set<int> set2);
        std::set<int> takeIntrsec(std::set<int>* setArray);
        bitArray query (std::string query_key, int len);
        void createMetaRambo(vector<string> setIDs, int K, bool verbose);
        void serializeRAMBO(std::string dir);
        void deserializeRAMBO(std::string dir);
        // void insertion2 (std::vector<std::string> alllines);
	bitArray queryseq (std::string query_key, int len);
	void insertionRare (std::string setID, std::vector<std::string> keys);
        // void insertionwithRead (std::string setID, std::string filenameSet);

        int R;
        int B;
        int n;
        float p;
        int range;
        int k;
        float FPR;
        IDLBloomFilter** Rambo_array;
        Hasher **hasher;
        std::vector<int>* metaRambo;
        std::vector<std::vector<unsigned int>> idx_and_r_to_b;
        unsigned int idx_file_to_b[Ki][Ri];
};

#endif
