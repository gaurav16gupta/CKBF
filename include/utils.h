#ifndef _utils_
#define _utils_
#include <vector>
#include <set>
#include <string>
#include <map>
#include <cassert>

std::vector<uint> myhash(std::string key, int len, int k, int range);
std::vector<uint> myhash2( std::string &&key, int len, int k, int range);
std::vector<uint> xxhash32( std::string &&key, int len, int k, int range);
uint *rand_nums(int len, int range);
std::vector<uint> tabulation_hash(std::string &&key, int k, uint rand_nums[]);
std::vector<uint> myhashCacheLoc(std::string key, int len, int k, int range);
std::vector<uint> myhashCheap(std::string key, int len, int k, int range);

std::vector <std::string> getFastqdata(std::string filenameSet);
std::vector<std::string> readlines( std::string path, int num);
void murmurhash( const char *  key, int len, int k, uint range, uint * hashvals);
// std::vector<std::string> getsets( std::string path);
// std::vector<std::string> line2array( std::string line, char d);
// void writeRAMBOresults(std::string path, int rows, int cols, float* values);
// std::vector<int>  arrayunion(std::vector<int> &v1, std::vector<int> &v2);
// std::vector<int> arrayintersection(std::vector<int> &v1, std::vector<int> &v2);
// std::set<int> takeunion(std::set<int> set1, std::set<int> set2);
// std::vector <std::string> getctxdata(std::string filenameSet);

// std::vector<std::string> getRandomTestKeys(int keysize, int n);
// std::map<std::string, std::vector<int>> makeInvIndex(int n, std::vector<std::string> foldernames);
// std::vector<std::string> getkmers(std::string query_key, int kmersize);

#endif
