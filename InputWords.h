//
// Created by henrique on 04/05/2020.
//

#ifndef GRAMMARINDUCTION_INPUTWORDS_H
#define GRAMMARINDUCTION_INPUTWORDS_H


#include <vector>
#include "Symbol.h"
#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <unordered_map>

namespace fs = std::experimental::filesystem;

using namespace std;

class InputWords {
public:
    vector<vector<Symbol>> inputWords;
    vector<vector<Symbol>> testWords;
    int actualShare;
private:
    bool timed;
    std::unordered_map<string, int> chordCounts;
    std::unordered_map<string, Symbol> chordMap;
    int nTerminals;
    int nTestShares;

public:
    InputWords(bool timed, int nTerminals);
    void readWords();
    void convertFiletoWord (fs::path path, int minSize);
    void transposeTo(string actualTone, string targetTone, vector<Symbol> & word);
    vector<Symbol> generateTerminals();
    void printWordSizes();
    void selectTrainingWords(int nSharesOrAmount, bool byShare);
    void countChords();
    void truncChordWords();
    void checkEmptyString(string path);
    bool nextShareTrainingWords();
};



#endif //GRAMMARINDUCTION_INPUTWORDS_H
