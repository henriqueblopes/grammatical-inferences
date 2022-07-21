#include "InputWords.h"
#include <algorithm>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <cfloat>
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
void load_spice_file(string filename, vector<Symbol::Symbol> &terminals, vector<vector<Symbol::Symbol>> &words_to_infer, int index, int n_symbol_max);
void chord_to_char(vector<Symbol::Symbol> &terminals, vector<vector<Symbol::Symbol>> &words_to_infer, vector<Symbol::Symbol> &char_symbols, vector<vector<Symbol::Symbol>> &char_words);
vector<double> load_pautomac_solution_file(string filename, int index);
vector<vector<Symbol::Symbol>> generate_palindromes(int max_length);
vector<vector<Symbol::Symbol>> generate_mod_a_eq_mod_b(int max_length);
vector<vector<Symbol::Symbol>> generate_expression_language(int max_length);


int main(int argc, char** argv) {
    if (argc < 9)
        exit(1);
    srand(time(0));



    //For ICAART
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

    /*int n_terminals = stoi((argv[1]));
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
        timed = false;*/



    std::vector<Symbol::Symbol> expression_terms = {Symbol::Symbol("(", 0, true, false), Symbol::Symbol(")", 1, true,false), Symbol::Symbol("a", 2, true,false), Symbol::Symbol("+", 3, true,false)};
    std::vector<Symbol::Symbol> terms = {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true,false) /*};*/, Symbol::Symbol("2", 2, true,false)};
    vector<vector<Symbol::Symbol>> words, chord_words;
    /*words = {{Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                    {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true,false),  Symbol::Symbol("1", 1, true,false),  Symbol::Symbol("1", 1,true,false)},
                                                                            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1,true,false)}};*/

     //READ MUSICAL DATABASE
    InputWords iw = InputWords(false, 50);
    iw.read_words();
    //iw.input_words.erase(iw.input_words.begin(), iw.input_words.begin()+995);
    iw.iterate_chords();
    iw.change_words_to_reducted_chords();
    vector<Symbol::Symbol> chord_terms = iw.generate_terminals(iw.reducted_chord_counts);
    chord_words = iw.input_words;
    chord_to_char(chord_terms, chord_words, terms, words);


    /* CREATE |a| = |b|
    Grammar::Grammar g_mod_a_mod_b = Grammar::Grammar(terms, 1, mod_a_eq_mod_b, g_mod_a_mod_b.pfa, make_pair(0, 0));
    vector<Symbol::Symbol> lhs;
    lhs.push_back(g_mod_a_mod_b.non_terminals[0]);
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right;
    vector<Symbol::Symbol> rhs;
    rhs.push_back(terms[0]); rhs.push_back(g_mod_a_mod_b.non_terminals[0]); rhs.push_back(terms[1]); rhs.push_back(g_mod_a_mod_b.non_terminals[0]);
    right.push_back(make_pair(rhs, make_pair(1, 0))); rhs.clear();
    rhs.push_back(terms[1]); rhs.push_back(g_mod_a_mod_b.non_terminals[0]); rhs.push_back(terms[0]); rhs.push_back(g_mod_a_mod_b.non_terminals[0]);
    right.push_back(make_pair(rhs, make_pair(1, 0))); rhs.clear();
    rhs.push_back(Symbol::Symbol("", -1, true, false));
    right.push_back(make_pair(rhs, make_pair(1, 0))); rhs.clear();
    Rule::Rule raux = Rule::Rule(lhs, right);
    g_mod_a_mod_b.rules.clear(); g_mod_a_mod_b.rules.push_back(raux);
    words = g_mod_a_mod_b.generate_max_size_words_from_rules(16);*/


    //For MLJ
    //argv[0] = train method 0 alergia, 1 PInference, 2 MH (int)
    //argv[1] = alhpa for alergia and PInference (double)
    //argv[2] = ratio for PInference
    //argv[3] = n_Nonterminals for MH
    //argv[4] = time limite (for all?)
    //argv[5] = n_iterations for MH
    //argv[6] = database (6, 14, 16, 18)
    //argv[7] = max word lenght (20 for 6 and 14 spice, 99 for 16, ? for 18)
    // pumpns achados em spice: 6, 14,
    int training_method = stoi((argv[1]));
    double alpha = stod(argv[2]);
    double p_ratio = stod(argv[3]);
    int n_non_terminals = stoi(argv[4]);
    double time_limit = stod(argv[5]);
    int iterations = stoi((argv[6]));
    int index_p_file = stoi((argv[7]));
    int max_word_lenght = stoi((argv[8]));



    //load_pautomac_file(".pautomac.train", terms, words, index_p_file);
    //load_spice_file(".spice.train", terms, words, index_p_file, max_word_lenght);
    double avg_score = 0.0;
    training_method = 1;
    for (int i = 0; i < 1; i ++) {
        cout << "it: " << i+1 << ": ";



        //std::random_shuffle ( words.begin(), words.end());
        vector<vector<Symbol::Symbol>> test_words;
        vector<vector<Symbol::Symbol>> train_words;
        test_words.insert(test_words.end(), words.begin(), words.begin()+words.size()/10);
        train_words.insert(train_words.end(), words.begin()+words.size()/10,  words.end());
        train_words = words;
        Grammar::Grammar g = Grammar::Grammar(terms, n_non_terminals, train_words, g.pcsg, make_pair(0, 0));
        /*for (auto w: g.words)
            cout << g.convert_vector_to_string(w) << endl;
        exit(0);*/
        if (training_method == 0) {
            g.g_tp = g.pfa;
        } else if (training_method == 1)
            g.g_tp = g.pcfg;
        else if (training_method == 2)
            g.g_tp = g.pcsg;
        else if (training_method == 3)
            g.g_tp = g.n_gram;
        else
            exit(-1);

        vector<double> sol_pal;
        for (auto w: test_words) {
            double count = 0;
            for (auto w2: test_words) {
                if (g.equal_word(w, w2)) {
                    count += 1.0;
                }
            }
            sol_pal.push_back(1.0 * count/(1.0 * test_words.size()));
        }

        if (training_method == 0)
            g.train(g.pfa_alergia, iterations, alpha, p_ratio, time_limit);
        else if (training_method == 1)
            g.train(g.pcfg_pumping_inference, iterations, alpha, p_ratio, time_limit);
        else if (training_method == 2)
            g.train(g.pcsg_metropolis_hastings, iterations, alpha, p_ratio, time_limit);

        long double exp = 0.0;
        long double expI = 0.0;
        long double pcx = 0.0;
        for (int i = 0 ; i < test_words.size(); i++) {
            if (training_method == 3)
                pcx = 1/ (1.0 * pow(terms.size()+1, test_words[i].size()));
            else if (training_method == 2)
                pcx = g.find_word_probabilities_from_pcfg_inside_table(test_words[i]);
            else
                pcx = g.find_word_probabilities(test_words[i]);

            //cout << "word "<< i << " prob: "<< pcx <<" - probSol:  " << sol_pal[i] << endl;
            if (pcx == 0.0)
                pcx = 1/ (1.0 * train_words.size() * (1.0 * pow(terms.size()+1, test_words[i].size())));

            exp += sol_pal[i] * log2(pcx);
            expI += sol_pal[i] * log2(sol_pal[i]);
        }
        avg_score += pow(2, -exp);
        cout << " Score: " << pow(2, -exp) << " IdealScore: " << pow(2, -expI)<< endl;
        g.print_grammar();
    }
    cout << "AVG Score: " << avg_score/30;

     //save_grammar(g, "grammar.txt");
    //Carregar Gramatica e gerar mapas
    //load_grammar(g, "grammar.txt");

    /* vector<pair<int, double>> v = g.find_prefix_ranking_probabilities(prefix);
    for (auto p: v)
        cout << " terminal: " << p.first << " prob " << p.second << endl;*/
    exit(0);





    /*pair<int,int> context_size;
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
        Grammar::Grammar g = Grammar::Grammar(terms, 3, words, g.pfa, make_pair(0, 0));
        g.train(g.pfa_collapsed_gibbs_sample, iterations);

        g.print_grammar();

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        pair<double,double> p1;
        cout  << " Calcutating perplexity..." << endl;
        *//*if (training_method == 3 || training_method == 1)
            p1 = g.perplexity(iw.test_words);
        else if (training_method == 4 || training_method == 2)
            p1 = g.perplexity_kl(iw.test_words);*//*

        cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
        cout << "Share: " << iw.actual_share << " - Final Perplexity: " << p1.first << " - Final NPerplexity: " << p1.second << endl;
        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        cout << "elapsed time with perplexity time: " << elapsed_seconds.count() << "s" << endl  << endl;
    } while (iw.next_share_training_words());*/



}

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

void load_spice_file(string filename, vector<Symbol::Symbol> &terminals, vector<vector<Symbol::Symbol>> &words_to_infer, int index, int n_symbol_max) {
    terminals.clear();
    words_to_infer.clear();
    fs::path p = fs::current_path().parent_path();
    fs::path path = (p /= "resources/SPiCe_Offline/train/" + to_string(index) + filename);
    std::ifstream ifs (path, std::ifstream::in);
    string line;
    getline(ifs,line);
    size_t tokenPos = line.find(" ");
    int n_terminals = stoi(line.substr(tokenPos,string::npos));
    if (n_terminals >= 71) {
        cout << "Error. Too many terminals. Max allowed is 71 terminals" << endl;
        return;
    }
    for (int i = 0; i < n_terminals; i++) {
        char symbol_name = 65;
        if (i < 10000)
            terminals.push_back(Symbol::Symbol( to_string(i), i, true, false));
        /*else {
            symbol_name += i;
            terminals.push_back(Symbol::Symbol( string(1, symbol_name), i, true, false));
        }*/
    }
    do  {
        getline(ifs,line);
        size_t tokenPos = line.find(" ");
        int n_symbol = 0;
        if (tokenPos != string::npos)
            n_symbol = stoi(line.substr(0,tokenPos+1));
        vector<Symbol::Symbol> word;

        if (n_symbol <=  n_symbol_max && n_symbol > 0 ) {
            line = line.substr(tokenPos+1, string::npos);
            while (tokenPos != string::npos) {

                tokenPos = line.substr(0, string::npos).find(" ");
                int symbol = stoi(line.substr(0,tokenPos));
                word.push_back(terminals[symbol]);
                line = line.substr(tokenPos+1, string::npos);
            }
        //if (!word.empty())
            words_to_infer.push_back(word);
        }



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
                //pal.insert(pal.end(), bpal.begin(), bpal.end());
                //pal.insert(pal.end(), bpal.begin(), bpal.end());
                pal.insert(pal.end(), aux.begin(), aux.end());
                pal.insert(pal.end(), bpal.begin(), bpal.end());
                //pal.insert(pal.end(), bpal.begin(), bpal.end());
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
            pal.insert(pal.end(), moda_0.begin(), moda_0.end());
            pal.insert(pal.end(), modb_1.begin(), modb_1.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            next_moda_modb_aux.push_back(pal);
            pal.clear();
            pal.insert(pal.end(), moda_0.begin(), moda_0.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            pal.insert(pal.end(), modb_1.begin(), modb_1.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            next_moda_modb_aux.push_back(pal);
            pal.clear();
            pal.insert(pal.end(), modb_1.begin(), modb_1.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            pal.insert(pal.end(), moda_0.begin(), moda_0.end());
            next_moda_modb_aux.push_back(pal);
            pal.clear();
            pal.insert(pal.end(), modb_1.begin(), modb_1.end());
            pal.insert(pal.end(), moda_0.begin(), moda_0.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            next_moda_modb_aux.push_back(pal);
            pal.clear();
            pal.insert(pal.end(), modb_1.begin(), modb_1.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
            pal.insert(pal.end(), moda_0.begin(), moda_0.end());
            pal.insert(pal.end(), aux.begin(), aux.end());
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
void chord_to_char(vector<Symbol::Symbol> &terminals, vector<vector<Symbol::Symbol>> &words_to_infer, vector<Symbol::Symbol> &char_symbols, vector<vector<Symbol::Symbol>> &char_words) {
    char_symbols.clear();
    for (int i = 0; i < terminals.size(); i++ ) {

        char c = 40 + i;
        char_symbols.push_back(Symbol::Symbol(std::string(1,c), terminals[i].id, terminals[i].terminal, terminals[i].context));
    }
    char_words = words_to_infer;
    for (auto & w: char_words)
        for (auto & t: w)
            t = char_symbols[t.id];
}
