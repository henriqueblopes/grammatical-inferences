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
    std::unordered_map<string, string> chordMapReduction;
    std::unordered_map<string, int> reductedChordCounts;
    int nTerminals;
    int nTestShares;

public:
    InputWords(bool timed, int nTerminals);
    void readWords();
    void convertFiletoWord (const fs::path& path, unsigned long minSize);
    static void transposeTo(const string& actualTone, const string& targetTone, vector<Symbol> & word);
    vector<Symbol> generateTerminals();
    void printWordSizes();
    void selectTrainingWords(int nSharesOrAmount, bool byShare);
    void countChords();
    void truncChordWords();
    void checkEmptyString(const string& path);
    bool nextShareTrainingWords();
    void iterateChords();
    static string buildChordVector(const string& mode, const string& addition, const string& bass);
private:
    static string addChordNote(string chord, string note);

    void countAndShowReductedChords();
};



#endif //GRAMMARINDUCTION_INPUTWORDS_H
