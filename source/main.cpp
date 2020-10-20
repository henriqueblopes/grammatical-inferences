#include <iostream>
/*
#include "Symbol.h"
#include "Grammar.h"
*/
#include <vector>
#include "InputWords.h"
#include "gInfer.h"



using namespace std;


int main(int argc, char** argv) {
    if (argc < 11)
        exit(1);
    //argv[1] = nTerminals
    //argv[2] = nSharesOrAmount
    //argv[3] = byShare
    //argv[4] = nInputForTraining
    //argv[5] = contextleftSize
    //argv[6] = contextRIghtSize
    //argv[7] = iterations
    //argv[8] = trainingMethod
    //argv[9] = normalizedPerplexity
    //argv[10] = timedInt
    //argv[11] = nNonterminals
    int nTerminals = stoi((argv[1]));
    int nSharesOrAmount = stoi((argv[2]));
    unsigned long nInputForTraining = stoi((argv[4]));
    int contextLeftSize = stoi((argv[5]));
    int contextRightSize = stoi((argv[6]));
    int iterations = stoi((argv[7]));
    int trainingMethod = stoi((argv[8]));
    int nNonterminals = stoi(argv[11]);
    bool byShare = true;
    if(!stoi((argv[3])))
        byShare = false;
    bool normalizedPerplexity = true;
    if (!stoi((argv[9])))
        normalizedPerplexity = false;
    bool timed = true;
    if(!stoi(argv[10]))
        timed = false;

    InputWords iw = InputWords(timed, nTerminals);
    iw.readWords();
    iw.iterateChords();

    std::vector<Symbol> terms = {Symbol("0", 0, true, false), Symbol("1", 1, true,false), Symbol("2", 2, true, false)};
    vector<vector<Symbol>> words = {{Symbol("1", 1, true, false), Symbol("1", 1, true, false), Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                    {Symbol("0", 0, true, false), Symbol("0", 0, true,false),  Symbol("1", 1, true,false),  Symbol("1", 1,true,false)},
                                    {Symbol("0", 0, true, false), Symbol("1", 1,true,false)}};

    vector<vector<Symbol>> words2 = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 0, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},

                                     {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},

                                     {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                     {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)}};

    //Grammar g = Grammar(terms, 3, words2, 4);
    //Grammar g = Grammar(terms, words2, 3);
    //g.printGrammar();
    //g.trainNGram();
    //g.printGrammar();


    //g.insideOutside(20);
    //g.printGrammar();

    vector<Symbol> terminals = iw.generateTerminals();
    iw.selectTrainingWords(nSharesOrAmount, byShare);
    if(iw.inputWords.size() < nInputForTraining)
        nInputForTraining = iw.inputWords.size();

    pair<int,int> contextSize;
    contextSize = make_pair(contextLeftSize,contextRightSize);

    cout << "nTerminals: "<< nTerminals << " " << "nSharesOrAmount: " << nSharesOrAmount << " byShare: " << byShare << " nInputForTraining: " <<
            nInputForTraining << " contextLeftSize: "  << contextLeftSize << endl<< " contextRightSize: " <<
            contextRightSize << " iterations: " << iterations << " trainingMethod: " << trainingMethod << endl << " normalizedPerplexity: " <<
            normalizedPerplexity << " timed: " << timed << " nNonterminals: " << nNonterminals << endl;

    do {
        auto start = std::chrono::system_clock::now();
        vector<vector<Symbol>> iWordsLim;
        iWordsLim.reserve(nInputForTraining);
        for (unsigned long i = 0; i < nInputForTraining; i++)
            iWordsLim.push_back(iw.inputWords[i]);

        Grammar g = Grammar(terminals, contextSize, nNonterminals, iWordsLim);
        //Grammar g = Grammar(terminals, contextSize, nNonterminals, iWordsLim, trainingMethod);

        g.train(trainingMethod,iterations);

        g.printGrammar();

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        pair<double,double> p1;
        cout  << " Calcutating perplexity..." << endl;
        if (trainingMethod == 3 || trainingMethod == 1)
            p1 = g.perplexity(iw.testWords);
        else if (trainingMethod == 4 || trainingMethod == 2)
            p1 = g.perplexityKL(iw.testWords);

        cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
        cout <<"Share: " << iw.actualShare<<" - Final Perplexity: " <<p1.first <<" - Final NPerplexity: " <<p1.second << endl;
        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        cout << "elapsed time with perplexity time: " << elapsed_seconds.count() << "s" << endl  << endl;
    } while (iw.nextShareTrainingWords());



}

/*std::string sentences[] = {
            "01",
            "000111",
            "00001111",
            "0011",
            "0000011111",
            "0101",
            "001101",
            "001011",
            "00011011",
            "010011",
            "001101",
            "010101",
            "00100111",
            "010011",
            "00110011",
            "000111",
            "0011",
            "0011001101",
            "01001101",
            "0101010101",
            "0001110011",
            "0011000111",
            "00101010101010101010010101010101010101010101"
    };*/

