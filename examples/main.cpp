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
void load_pautomac_file(string filename, vector<Symbol::Symbol> & terminals, vector<vector<Symbol::Symbol>> &words, int index);
void load_spice_file(string filename, vector<Symbol::Symbol> & terminals, vector<vector<Symbol::Symbol>> &words, int index);
vector<double> load_pautomac_solution_file(string filename, int index);
vector<vector<Symbol::Symbol>> generate_palindromes(int max_length);
vector<vector<Symbol::Symbol>> generate_mod_a_eq_mod_b(int max_length);
vector<vector<Symbol::Symbol>> generate_expression_language(int max_length);


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

    std::vector<Symbol::Symbol> expression_terms = {Symbol::Symbol("(", 0, true, false), Symbol::Symbol(")", 1, true,false), Symbol::Symbol("a", 2, true,false), Symbol::Symbol("+", 3, true,false)};
    std::vector<Symbol::Symbol> terms = {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true,false) /*};*/, Symbol::Symbol("2", 2, true,false)};
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

    vector<vector<Symbol::Symbol>> mod_a_eq_mod_b = generate_mod_a_eq_mod_b(15);

    vector<vector<Symbol::Symbol>> anbn = {
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)},
            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("0", 0, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false),Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false)}
    };
    vector<vector<Symbol::Symbol>> palindromes = generate_palindromes(30);
    vector<vector<Symbol::Symbol>> reverse_palindromes;
    for (auto p: palindromes) {
        vector<Symbol::Symbol> rp;
        for (int i = p.size()-1; i >=0; i--) {
            rp.push_back(p[i]);
        }
        reverse_palindromes.push_back(rp);
    }
    vector<vector<Symbol::Symbol>> expressions = generate_expression_language(25);



    /*Grammar::Grammar g = Grammar::Grammar(terms, 3, words2, g.pfa, make_pair(0, 0));
    g.train(g.pfa_alergia, 10);*/
    //iw.input_words.erase(iw.input_words.begin()+5790,iw.input_words.end());
    //iw.input_words.erase(iw.input_words.begin(),iw.input_words.begin()+5010);



    /*iw.iterate_chords();
    iw.change_words_to_reducted_chords();
    vector<Symbol::Symbol> music_terminals = iw.generate_terminals(iw.reducted_chord_counts);
    Grammar::Grammar g = Grammar::Grammar(music_terminals, 3, iw.input_words, g.pcfg, make_pair(0, 0));*/

    int index_p_file = 3;

    load_pautomac_file(".pautomac.train", terms, words2, index_p_file);
    //load_spice_file(".spice.train", terms, words2, index_p_file);
    //words2.erase(words2.begin()+1000, words2.end());

    Grammar::Grammar g = Grammar::Grammar(terms, 3, words2, g.pfa, make_pair(0, 0));
    /*g.train(Grammar::Grammar::pfa_alergia, 50);
    //g.print_grammar();
    exit(0);*/
    /*for (auto nt: g.non_terminals)
        for (auto nt2: g.non_terminals)
            if (g.compatible_alergia(nt, nt2, 0.95, g.rules))
                cout << nt.name << " and " << nt2.name << " compatible: " << g.compatible_alergia(nt, nt2, 0.95, g.rules) << endl;*/
    //g.words.erase(g.words.begin()+1000, g.words.end());

    /*for (auto e: g.words)
        cout << g.convert_vector_to_string(e) << endl;
    cout << "size: " <<g.words.size() << endl;*/


                // PUMPING INFERENCE BELOW //
    unordered_map<std::string, std::vector<int>> map_pump_to_word;
    unordered_map<std::string, int> map;
    unordered_map<int, vector<string>> node_pump_map;

/*
    vector<vector<Symbol::Symbol>> test_words;
    vector<Symbol::Symbol> aux;
    load_pautomac_file(".pautomac.test", aux, test_words, index_p_file);
    vector<double> sol = load_pautomac_solution_file(".pautomac_solution.txt", index_p_file);*/
    /*double e = 0.0;
    double psol;
    for (int i = 0 ; i < test_words.size(); i++)
        e += sol[i] * log2(sol[i]);
    cout << "perplexitiyy sol: " << pow(2,-e);
    exit(0);*/
    //inferir para cnf e salvar
    /*g.prob_sequitur();
    g.convert_to_cnf();
    g.remove_unused_rules();
*/
    //exit(0);
    //save_grammar(g, "grammar.txt");

    //verifica se gramática gera base de dados
    /*for (int i = 0; i < g.words.size(); i++) {
        std::vector<Symbol::Symbol> w = g.yield_string(g.rules[0].right[i].first);
        if(g.equal_word(w,g.words[776]))
            cout << "word " << i << " OK!" << endl;
        else
            cout << "word " << i << " Fail!" << endl;
    }*/

    //Carregar Gramatica e gerar mapas
    //load_grammar(g, "grammar.txt");
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    /*cout <<"Starting to find Pumps by Slice at " << std::ctime(&end_time) << endl;
    for (int i = 0; i < g.words.size(); i++) {
        if ((i+1) % (g.words.size()/5) == 0 )
            cout << "   word: " << i + 1 << ". " << (100.0*(i+1)/(g.words.size()*1.0)) << "%" << endl;
        if (g.words[i].size() <= 999999) {
            //g.count_pumping_str(map, g.words[i], 0, map_pump_to_word, g.words[i]);]
            g.count_pumping_str_by_slice(map, g.words[i], map_pump_to_word, i, 0.0);
        }
    }

    for (auto i = map.begin(), j = map.end(); i!= j; ) {
        if (i->second < 2) {
            map_pump_to_word.erase(i->first);
            i = map.erase(i);
        } else
            ++i;
    }
    g.eliminate_covered_pumpings(map_pump_to_word, map);




    vector<pair<string,int>> ordered;
    for (auto w : map)
        ordered.push_back(w);
    std::sort(ordered.begin(), ordered.end(), sortbysecond);
    for(auto it = ordered.begin(); it != ordered.end(); ++it)
        std::cout << it->first << ": " << it->second << endl;*/
    cout << "size: " << map.size() << " size database: " << g.words.size() << endl;
    end = std::chrono::system_clock::now();
    end_time = std::chrono::system_clock::to_time_t(end);
    cout << "Starting to Infer at " << std::ctime(&end_time) << endl;
    //g2.fpta_pumping_inference(map, map_pump_to_word, 1, 1, 3, 0);
    g.super_duper_pumping_inference();
    end = std::chrono::system_clock::now();
    end_time = std::chrono::system_clock::to_time_t(end);
    cout << "Finished  to Infer at " << std::ctime(&end_time) << endl;
    //g.print_grammar();
    //cout << "Finished  to Infer at " << std::ctime(&end_time) << endl;
    /*for (int i = 0; i < g.words.size(); i++)
        cout << "Word " << i << ": " << g.convert_vector_to_string(g.words[i]) << endl;*/
    /*vector<Symbol::Symbol> prefix = g.words[18];
    prefix.erase(prefix.begin()+4, prefix.end());
    cout << "Prefix " << 20 << ": " << g.convert_vector_to_string(prefix) << endl;
    vector<pair<int, double>> v = g.find_prefix_ranking_probabilities(prefix);
    for (auto p: v)
        cout << " terminal: " << p.first << " prob " << p.second << endl;*/
    exit(0);

    /*for (auto s: ordered)
        cout << s.first << ": " << s.second << endl;*/

    /*save_map(map, "map.txt");
    save_map_pump_to_word(map_pump_to_word, "map_pump_to_word.txt");*/


    //Pumping Inference
    /*load_grammar(g, "grammar.txt");
    load_map_pump_to_word(map_pump_to_word, "map_pump_to_word.txt");
    load_map(map, "map.txt");*/
    //vector<pair<string,int>> ordered;
    //g.print_grammar();
    std::cout << "Before pump: Total Rules: " << g.rules.size() << " Total initial rules: " << g.rules[0].right.size() << endl;
    g.pumping_inference(map, map_pump_to_word);
    std::cout << "After pump: Total Rules: " << g.rules.size() << " Total initial rules: " << g.rules[0].right.size() << endl;
    int count = 0;
    for (const auto& r: g.rules[0].right) {
        if (r.second.first > 0.0)
            count ++;
    }
    cout << "Righties: " << count << endl;
    //g.print_grammar();

   /* long double exp = 0.0;
    long double pcx = 0.0;
    for (int i = 0 ; i < test_words.size(); i++) {
        pcx = g.calculate_parse_tree_prob_top_down(test_words[i]);
        cout << "word "<< i << " prob: "<< pcx <<" - probSol:  " << sol[i] << endl;
        exp += sol[i] * log2(pcx);
    }
    cout << "Score: " << pow(2, -exp) << endl;*/
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


    exit(0);
    //g.print_grammar();
    //g.train_n_gram();
    //g.print_grammar();
    //g.inside_outside(20);
    //g.print_grammar();


    vector<Symbol::Symbol> terminals = iw.generate_terminals(std::unordered_map<string, int>());
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
void load_pautomac_file(string filename, vector<Symbol::Symbol> & terminals, vector<vector<Symbol::Symbol>> &words_to_infer, int index) {
    terminals.clear();
    words_to_infer.clear();
    fs::path p = fs::current_path().parent_path();
    fs::path path = (p /= "resources/Pautomac/" + to_string(index) + filename);
    std::ifstream ifs (path, std::ifstream::in);
    string line;
    getline(ifs,line);
    size_t tokenPos = line.find(" ");
    int n_terminals = stoi(line.substr(tokenPos,string::npos));
    for (int i = 0; i < n_terminals; i++)
        terminals.push_back(Symbol::Symbol(to_string(i), i, true, false));
    do  {
        getline(ifs,line);
        size_t tokenPos = line.find(" ");
        line = line.substr(tokenPos+1, string::npos);
        vector<Symbol::Symbol> word;
        while (tokenPos != string::npos) {

            tokenPos = line.substr(0, string::npos).find(" ");
            int symbol = stoi(line.substr(0,tokenPos));
            word.push_back(terminals[symbol]);
            line = line.substr(tokenPos+1, string::npos);
        }
        if (!word.empty())
            words_to_infer.push_back(word);


    } while (ifs.good());
}

void load_spice_file(string filename, vector<Symbol::Symbol> & terminals, vector<vector<Symbol::Symbol>> &words_to_infer, int index) {
    terminals.clear();
    words_to_infer.clear();
    fs::path p = fs::current_path().parent_path();
    fs::path path = (p /= "resources/SPiCe_Offline/train/" + to_string(index) + filename);
    std::ifstream ifs (path, std::ifstream::in);
    string line;
    getline(ifs,line);
    size_t tokenPos = line.find(" ");
    int n_terminals = stoi(line.substr(tokenPos,string::npos));
    if (n_terminals >= 35) {
        cout << "Error. Too many terminals. Max allowed is 35 terminals" << endl;
        return;
    }
    for (int i = 0; i < n_terminals; i++) {
        char symbol_name = 87;
        if (i < 10)
            terminals.push_back(Symbol::Symbol( to_string(i), i, true, false));
        else {
            symbol_name += i;
            terminals.push_back(Symbol::Symbol( string(1, symbol_name), i, true, false));
        }
    }
    do  {
        getline(ifs,line);
        size_t tokenPos = line.find(" ");
        line = line.substr(tokenPos+1, string::npos);
        vector<Symbol::Symbol> word;
        while (tokenPos != string::npos) {

            tokenPos = line.substr(0, string::npos).find(" ");
            int symbol = stoi(line.substr(0,tokenPos));
            word.push_back(terminals[symbol]);
            line = line.substr(tokenPos+1, string::npos);
        }
        //if (!word.empty())
            words_to_infer.push_back(word);


    } while (ifs.good());
    words_to_infer.erase(words_to_infer.end()-1, words_to_infer.end());
}

vector<double> load_pautomac_solution_file(string filename, int index) {
    fs::path p = fs::current_path().parent_path();
    fs::path path = (p /= "resources/Pautomac/" + to_string(index) + filename);
    std::ifstream ifs (path, std::ifstream::in);
    string line;
    getline(ifs,line);
    vector<double> solutions;
    while (ifs.good()) {
        getline(ifs,line);
        if (!line.empty())
            solutions.push_back(stod(line));
    }
    return solutions;

}

vector<vector<Symbol::Symbol>> generate_palindromes(int max_length) {
    vector<Symbol::Symbol> pal_0 = {Symbol::Symbol("0", 0, true, false)};
    vector<Symbol::Symbol> pal_1 = {Symbol::Symbol("1", 1, true, false)};
    vector<Symbol::Symbol> pal_2 = {Symbol::Symbol("2", 2, true, false)};
    vector<Symbol::Symbol> pal_e = {};
    vector<vector<Symbol::Symbol>> all_palindromes;

    //all_palindromes.push_back(pal_e);
    //all_palindromes.push_back(pal_0);
    all_palindromes.push_back(pal_1);

    all_palindromes[0].push_back(pal_1[0]);
    vector<vector<Symbol::Symbol>> basic_palindromes;
    basic_palindromes.push_back(pal_0);
    basic_palindromes.push_back(pal_1);
    vector<vector<Symbol::Symbol>> next_palindromes = all_palindromes;
    vector<vector<Symbol::Symbol>> next_palindromes_aux;
    while (all_palindromes[all_palindromes.size()-1].size() < max_length) {
        for (auto bpal: basic_palindromes) {
            for (auto pal : next_palindromes) {
                vector<Symbol::Symbol> aux = pal;
                pal.clear();
                pal.insert(pal.end(), bpal.begin(), bpal.end());
                pal.insert(pal.end(), bpal.begin(), bpal.end());
                //pal.insert(pal.end(), bpal.begin(), bpal.end());
                pal.insert(pal.end(), aux.begin(), aux.end());
                pal.insert(pal.end(), bpal.begin(), bpal.end());
                pal.insert(pal.end(), bpal.begin(), bpal.end());
                pal.insert(pal.end(), bpal.begin(), bpal.end());
                next_palindromes_aux.push_back(pal);
            }
        }
        all_palindromes.insert(all_palindromes.end(), next_palindromes_aux.begin(), next_palindromes_aux.end());
        next_palindromes = next_palindromes_aux;
        next_palindromes_aux.clear();
    }
    next_palindromes = {pal_0};
    all_palindromes.push_back(pal_0);
    while (all_palindromes[all_palindromes.size()-1].size() < max_length) {
        for (auto pal : next_palindromes) {
            pal.insert(pal.end(), pal_2.begin(), pal_2.end());
            next_palindromes_aux.push_back(pal);
        }
        all_palindromes.insert(all_palindromes.end(), next_palindromes_aux.begin(), next_palindromes_aux.end());
        next_palindromes = next_palindromes_aux;
        next_palindromes_aux.clear();
    }
    for (auto & a: all_palindromes) {
        a.push_back(pal_2[0]);
        //reverse(a.begin(), a.end());
    }
    return all_palindromes;
}

vector<vector<Symbol::Symbol>> generate_mod_a_eq_mod_b(int max_length) {
    vector<Symbol::Symbol> moda_0 = {Symbol::Symbol("0", 0, true, false)};
    vector<Symbol::Symbol> modb_1 = {Symbol::Symbol("1", 1, true, false)};
    vector<Symbol::Symbol> pal_e = {};
    vector<vector<Symbol::Symbol>> all_moda_modb;

    all_moda_modb.push_back(pal_e);
    //all_moda_modb.push_back(moda_0);
    //all_moda_modb.push_back(modb_1);
    vector<vector<Symbol::Symbol>> basic_moda_modb;
    basic_moda_modb.push_back(moda_0);
    basic_moda_modb.push_back(modb_1);
    vector<vector<Symbol::Symbol>> next_moda_modb = all_moda_modb;
    vector<vector<Symbol::Symbol>> next_moda_modb_aux;
    while (all_moda_modb[all_moda_modb.size()-1].size() < max_length) {
        for (auto pal : next_moda_modb) {
            vector<Symbol::Symbol> aux = pal;
            pal.clear();
            pal.insert(pal.end(), moda_0.begin(), moda_0.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            pal.insert(pal.end(), modb_1.begin(), modb_1.end());
            next_moda_modb_aux.push_back(pal);
            pal.clear();
            pal.insert(pal.end(), modb_1.begin(), modb_1.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            pal.insert(pal.end(), moda_0.begin(), moda_0.end());
            next_moda_modb_aux.push_back(pal);
        }
        all_moda_modb.insert(all_moda_modb.end(), next_moda_modb_aux.begin(), next_moda_modb_aux.end());
        next_moda_modb = next_moda_modb_aux;
        next_moda_modb_aux.clear();
    }
    return all_moda_modb;
}

vector<vector<Symbol::Symbol>> generate_expression_language(int max_length) {
    vector<Symbol::Symbol> bracket_o = {Symbol::Symbol("(", 0, true, false)};
    vector<Symbol::Symbol> bracket_c = {Symbol::Symbol(")", 1, true, false)};
    vector<Symbol::Symbol> a = {Symbol::Symbol("a", 2, true, false)};
    vector<Symbol::Symbol> plus = {Symbol::Symbol("+", 3, true, false)};
    vector<Symbol::Symbol> pal_e = {};
    vector<vector<Symbol::Symbol>> all_expressions;

    //all_moda_modb.push_back(pal_e);
    all_expressions.push_back(a);
    vector<vector<Symbol::Symbol>> basic_moda_modb;

    vector<vector<Symbol::Symbol>> next_expressions = all_expressions;
    vector<vector<Symbol::Symbol>> next_expressions_aux;
    while (all_expressions[all_expressions.size()-1].size() < max_length) {
        for (auto pal : next_expressions) {
            vector<Symbol::Symbol> aux = pal;
            pal.clear();
            pal.insert(pal.end(), bracket_o.begin(), bracket_o.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            pal.insert(pal.end(), bracket_c.begin(), bracket_c.end());
            next_expressions_aux.push_back(pal);

            /*for (auto pal2 : all_expressions) {
                pal.clear();
                pal.insert(pal.end(), aux.begin(), aux.end());
                pal.insert(pal.end(), plus.begin(), plus.end());
                pal.insert(pal.end(), pal2.begin(), pal2.end());
                next_expressions_aux.push_back(pal);
            }*/
        }
        all_expressions.insert(all_expressions.end(), next_expressions_aux.begin(), next_expressions_aux.end());
        next_expressions = next_expressions_aux;
        next_expressions_aux.clear();
    }
    return all_expressions;
}