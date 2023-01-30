#ifndef _UTILS_H_
#define _UTILS_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

vector<string> getFastqData(string fileName) {
    vector<string> lines;
    ifstream fastqFile(fileName);
    string line;
    uint64_t line_num = 0;
    while (getline(fastqFile, line)) {
        if (line_num % 4 == 1 && line.size() >= 31) {
            lines.push_back(line);
        }
        line_num += 1;
    }
    return lines;
}

vector<string> getQueryData(string fileName) {
    vector<string> lines;
    ifstream queryFile(fileName);
    string line;
    uint64_t line_num = 0;
    while (getline(queryFile, line)) {
        if (line.size() >= 31) {
            lines.push_back(line);
        }
        line_num += 1;
    }
    return lines;
}

void getQueryForArBF(string queryFileName, vector<string> &queries) {
    ifstream queryFile(queryFileName);
    string line;
    while (getline(queryFile, line) && line.size() > 0) {
        queries.push_back(line);
    }
}

void getQueryforArBF(string queryFilename, vector<string>& queries, vector<vector<uint32_t>>& GT, uint32_t N) {
    ifstream queryFile(queryFilename);
    string line;
    uint32_t line_num = 0;
    while (getline(queryFile, line)) {
        stringstream is;
        is<<line;
        string vals;
        vector <string> op;
        bool first =true;
        while( getline (is, vals, ',')){
            if (first) {queries[line_num] = vals; first=false;}
            else {
                if (stoi(vals)<=N) GT[line_num].push_back(stoi(vals));   
            }
        }    
        line_num += 1;
    }
}


// Vector reader
// template<typename T>
vector<uint64_t> read_vector(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<uint64_t> result;
    uint64_t value;
    while (file >> value) {
        result.push_back(value);
    }
    return result;
}

// Vector writer 
void write_vector(const std::string& filename, const std::vector<uint64_t>& vector) { 
    std::ofstream file(filename); 
    for (const auto& element : vector) { 
        file << element << " "; 
    } 
    file.close(); 
}

#endif
