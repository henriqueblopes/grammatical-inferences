#include <iostream>
#include <vector>
#include "InputWords.h"
#include <gInfer/Grammar.h>



using namespace std;


int main(int argc, char** argv) {
    if (argc < 11)
        exit(1);
    //argv[1] = n_terminals
    //argv[2] = n_shares_or_amount
    //argv[3] = by_share
    //argv[4] = n_input_for_training
    //argv[5] = contextleftSize
    //argv[6] = contextRIghtSize
    //argv[7] = iterations
    //argv[8] = training_method
    //argv[9] = normalized_perplexity
    //argv[10] = timedInt
    //argv[11] = n_non_terminals
    int n_terminals = stoi((argv[1]));
    int n_shares_or_amount = stoi((argv[2]));
    size_t n_input_for_training = stoi((argv[4]));
    int context_left_size = stoi((argv[5]));
    int context_right_size = stoi((argv[6]));
    int iterations = stoi((argv[7]));
    int training_method = stoi((argv[8]));
    int n_non_terminals = stoi(argv[11]);
    bool by_share = true;
    if(!stoi((argv[3])))
        by_share = false;
    bool normalized_perplexity = true;
    if (!stoi((argv[9])))
        normalized_perplexity = false;
    bool timed = true;
    if(!stoi(argv[10]))
        timed = false;

    InputWords iw = InputWords(timed, n_terminals);
    iw.read_words();
    iw.iterate_chords();

    std::vector<Symbol::Symbol> terms = {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true,false), Symbol::Symbol("2", 2, true, false)};
    vector<vector<Symbol::Symbol>> words = {{Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                    {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true,false),  Symbol::Symbol("1", 1, true,false),  Symbol::Symbol("1", 1,true,false)},
                                    {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1,true,false)}};

    vector<vector<Symbol::Symbol>> words2 = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 0, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},{Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false)},

                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},

                                     {Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false)},
                                     {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false)}};


    //g.print_grammar();
    //g.train_n_gram();
    //g.print_grammar();
    //g.inside_outside(20);
    //g.print_grammar();

    vector<Symbol::Symbol> terminals = iw.generate_terminals();
    iw.select_training_words(n_shares_or_amount, by_share);
    if(iw.input_words.size() < n_input_for_training)
        n_input_for_training = iw.input_words.size();

    pair<int,int> context_size;
    context_size = make_pair(context_left_size, context_right_size);

    cout << "n_terminals: " << n_terminals << " " << "n_shares_or_amount: " << n_shares_or_amount << " by_share: " << by_share << " n_input_for_training: " <<
         n_input_for_training << " context_left_size: " << context_left_size << endl << " context_right_size: " <<
         context_right_size << " iterations: " << iterations << " training_method: " << training_method << endl << " normalized_perplexity: " <<
         normalized_perplexity << " timed: " << timed << " n_non_terminals: " << n_non_terminals << endl;

    do {
        auto start = std::chrono::system_clock::now();
        vector<vector<Symbol::Symbol>> i_words_lim;
        i_words_lim.reserve(n_input_for_training);
        for (unsigned long i = 0; i < n_input_for_training; i++)
            i_words_lim.push_back(iw.input_words[i]);
        Grammar::Grammar g = Grammar::Grammar(terms, 3, words2, g.pfa, make_pair(0, 0));
        g.train(g.pfa_collapsed_gibbs_sample, iterations);

        g.print_grammar();

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        pair<double,double> p1;
        cout  << " Calcutating perplexity..." << endl;
        /*if (training_method == 3 || training_method == 1)
            p1 = g.perplexity(iw.test_words);
        else if (training_method == 4 || training_method == 2)
            p1 = g.perplexity_kl(iw.test_words);*/

        cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
        cout << "Share: " << iw.actual_share << " - Final Perplexity: " << p1.first << " - Final NPerplexity: " << p1.second << endl;
        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        cout << "elapsed time with perplexity time: " << elapsed_seconds.count() << "s" << endl  << endl;
    } while (iw.next_share_training_words());



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

