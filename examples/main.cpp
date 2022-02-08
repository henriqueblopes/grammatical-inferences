#include "InputWords.h"
#include <algorithm>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <gInfer/Grammar.h>
#include <iostream>
#include <map>
#include <vector>


using namespace std;

bool sortbysecond(const pair<string ,int> &a,
               const pair<string,int> &b)
{
    return (a.second > b.second);
}

void save_grammar(Grammar::Grammar &g, std::string filename);
void load_grammar(Grammar::Grammar &g, std::string filename);
void save_map_pump_to_word(std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::string filename);
void load_map_pump_to_word(std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::string filename);
void save_map(unordered_map<std::string, int> &map, std::string filename);
void load_map(unordered_map<std::string, int> &map, std::string filename);

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
    //iw.input_words.erase(iw.input_words.begin(), iw.input_words.begin()+995);
    //iw.iterate_chords();

    std::vector<Symbol::Symbol> terms = {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true,false)};
    vector<vector<Symbol::Symbol>> words = {{Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                    {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true,false),  Symbol::Symbol("1", 1, true,false),  Symbol::Symbol("1", 1,true,false)},
                                                                            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1,true,false)}};
    vector<vector<Symbol::Symbol>> lastWord = { {Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("2", 2, true, false),Symbol::Symbol("3", 3, true, false), Symbol::Symbol("1", 1, true, false),Symbol::Symbol("2", 2, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("2", 2, true, false),Symbol::Symbol("3", 3, true, false)}};


    lastWord = { {Symbol::Symbol("a",0,true,false),//a
                  Symbol::Symbol("b",1,true,false),//b
                  Symbol::Symbol("d",2,true,false),//d
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  //Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("d",2,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("s",3,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("b",1,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("b",1,true,false),
                  Symbol::Symbol("a",0,true,false),
                  Symbol::Symbol("b",1,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("g",5,true,false),
                  Symbol::Symbol("f",4,true,false),
                  Symbol::Symbol("g",5,true,false)}};





    vector<vector<Symbol::Symbol>> words2 = {/*{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                             {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                             {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                             {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                             {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},*/
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


    /*Grammar::Grammar g = Grammar::Grammar(terms, 3, words2, g.pfa, make_pair(0, 0));
    g.train(g.pfa_alergia, 10);*/
    //iw.input_words.erase(iw.input_words.begin()+5,iw.input_words.end());

    vector<Symbol::Symbol> music_terminals = iw.generate_terminals();
    Grammar::Grammar g = Grammar::Grammar(music_terminals, 3, iw.input_words, g.pcfg, make_pair(0, 0));
    unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> map_pump_to_word;
    unordered_map<std::string, int> map;

    //inferir para cnf e salvar
    /*g.prob_sequitur();
    g.convert_to_cnf();
    g.remove_unused_rules();*/
    //exit(0);
    //save_grammar(g, "grammar.txt");

    //verifica se gramática gera base de dados
    /*for (int i = 0; i < g.words.size(); i++) {
        std::vector<Symbol::Symbol> w = g.yield_string(g.rules[0].right[i].first);
        if(g.equal_word(w,g.words[i]))
            cout << "word " << i << " OK!" << endl;
        else
            cout << "word " << i << " Fail!" << endl;
    }*/

    //Carregar Gramatica e gerar mapas
    load_grammar(g, "grammar.txt");
    for (int i = 835; i < g.words.size(); i++) {
        cout << "word: " << i + 1 << endl;
        if (g.words[i].size() <= 999999) {
            //g.count_pumping_str(map, g.words[i], 0, map_pump_to_word, g.words[i]);
            g.count_pumping_str_by_slice(map, g.words[i], map_pump_to_word);
        }
    }
    /*save_map(map, "map.txt");
    save_map_pump_to_word(map_pump_to_word, "map_pump_to_word.txt");*/


    //Pumping Inference
    /*load_grammar(g, "grammar.txt");
    load_map_pump_to_word(map_pump_to_word, "map_pump_to_word.txt");
    load_map(map, "map.txt");*/
    //vector<pair<string,int>> ordered;
    //g.pumping_inference(map, map_pump_to_word);
    //save_grammar(g, "grammarPumped.txt");

    /*g.count_pumping_str(map, iw.input_words[21], 0);
    for (int i = 0; i < iw.input_words.size(); i++) {
        cout << "word: " << i+1 << endl;
        cout << "chords: "+ to_string(iw.input_words[i].size()) + " - ";
        for (auto a: iw.input_words[i])
            cout << a.name + " ";
        cout << endl;
        //30 é muito
        if (iw.input_words[i].size() < 24)
            g.count_pumping_str(map, iw.input_words[i], 0);
    }*/

    vector<pair<string,int>> ordered;
    for (auto w : map)
        ordered.push_back(w);
    std::sort(ordered.begin(), ordered.end(), sortbysecond);
    for(auto it = ordered.begin(); it != ordered.end(); ++it)
        std::cout << it->first << ": " << it->second << endl;

    /*for (auto s: ordered)
        cout << s.first << ": " << s.second << endl;*/
    cout << "size: " << map.size() << " size database: " << g.words.size() << endl;
    exit(0);
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

void save_grammar(Grammar::Grammar &g, std::string filename) {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << g;
}

void load_grammar(Grammar::Grammar &g, std::string filename) {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> g;

}

void save_map_pump_to_word(std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::string filename) {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << map_pump_to_word;
}

void load_map_pump_to_word(std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::string filename) {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> map_pump_to_word;
}

void save_map(unordered_map<std::string, int> &map, std::string filename) {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << map;
}

void load_map(unordered_map<std::string, int> &map, std::string filename) {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> map;
}