#ifndef _UTILS_H_
#define _UTILS_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

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

#endif