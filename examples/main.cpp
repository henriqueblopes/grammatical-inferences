#include "InputWords.h"
#include "midi_handler/MidiFile.h"
#include "midi_handler/Options.h"
#include <algorithm>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <cfloat>
#include <gInfer/Grammar.h>
#include <iostream>
#include <map>
#include <vector>
/*My SPiCe
 * 20 Billboard Songs
 * 21 Billboard Phrases
 * 22 Brown Tokens 10%
 * 23 Brown Tokens Complete
 * 24 Brown Part of Speech Complete (Acha PumpDet Duplo)
 * 25 is tail embedding
 * 26 Dyck Language 2_10 (Acha Pump duplo)
 * 27 Dyck Language 2_12 (Acha Pump duplo)
 * 28 Dyck Language 2_14 (Acha Pump
 * 29 Dyck Language 2_16
 * 30 Dyck Language 3_8
 * 31 Dyck Language 3_10 (Acha Pump)
 * 32 Dyck Language 3_12
 * 33 Dyck Language 4_8
 * 34 Dyck Language 4_10 (Acha Pump)
 * 35 Dyck Language 5_8 (Acha Pump
 * 36 Dyck Language 5_10
 * 37 Dyck Language 6_6
 * 38 Dyck Language 6_8
 * 39 Dyck Language 6_10
 * 40 Dyck Language 4_12
 * 41 Conll2003 POS full:14987 11:8751 21:11289
 * 42 Conll2003 POS with Valid full: 18450 11:10659 21: 13789
 * */


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
void write_string_to_midi(vector<Symbol::Symbol> string, std::string filename);
vector<double> load_pautomac_solution_file(string filename, int index);
vector<vector<Symbol::Symbol>> generate_palindromes(int max_length);
vector<vector<Symbol::Symbol>> generate_mod_a_eq_mod_b(int max_length);
vector<vector<Symbol::Symbol>> generate_expression_language(int max_length);
vector<vector<Symbol::Symbol>> generate_tail_embedding_language(int max_length);
vector<vector<Symbol::Symbol>> generate_dyck_n_embedding_language(int n, int max_length);
/*TODO
 * gerações Conll 10 20 30 size word - 2,3,5,7,10,15 NTs, 1000 iterações
 *
 * Fazer o Brown com POS
 * Experimento com spice 13 NLP (English spelling correction from Twitter Typos Corpus) (achou pumping v14 e v89)
 * Generate Dyck Languges ()[]{} length less than 55, 6230 words, test with 1000 words with length 20<l<70 (d6 have 15000 for Train,  2000 for test)
 * And this grammar (tail embedding english)
 *  S→NP VP
 *  VP→V1 | V2 NP
 *  NP→N | N CP
 *  CP→R VP
 *
 *  */



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
    /*InputWords iw = InputWords(false, 70);
    //iw.read_words_beatles(true);
    //iw.input_words.erase(iw.input_words.begin(), iw.input_words.begin()+995);
    iw.iterate_chords();
    iw.change_words_to_reducted_chords();
    vector<Symbol::Symbol> chord_terms = iw.generate_terminals(iw.reducted_chord_counts);
    chord_words = iw.input_words;
    cout << chord_terms.size() << endl;
    chord_to_char(chord_terms, chord_words, terms, words);*/



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
    //index_p_file = 23;

    //words = generate_dyck_n_embedding_language(index_p_file,max_word_lenght);
    /*InputWords iw = InputWords(false, INT64_MAX);
    iw.read_words_conll2003();
    iw.count_chords();
    vector<Symbol::Symbol> chord_terms = iw.generate_terminals(iw.chord_counts);
    chord_terms = iw.change_symbols_to_numbers(chord_terms);
    words = iw.input_words;
    ofstream myfile;
    myfile.open ("42.spice.train");
    myfile << words.size();
    myfile << " " << chord_terms.size() << endl;
    for (auto w: words) {
        myfile << w.size();
        for (auto s: w)
            myfile << " " << s.name;
        myfile << endl;
    }
    myfile.close();
    exit(1);*/
    /*
    Grammar::Grammar g = Grammar::Grammar(chord_terms, 3, chord_words, g.pcfg, make_pair(0, 0));
    g.train(g.pcfg_pumping_inference, iterations, alpha, p_ratio, time_limit);
    //load_grammar(g, "cfdg_grammar_full_beatles1.500000_0.200000_21.txt");
    //g.print_grammar();
    exit(0);*/


    /*for (int i = 0; i < 565; i++)
        cout << "0.1 NT0 --> " << to_string(i) << endl;
    exit(0);*/


    //load_pautomac_file(".pautomac.train", terms, words, index_p_file);
    load_spice_file(".spice.train", terms, words, index_p_file, max_word_lenght);
    
    //contar termiunais usados
    /*int count_t = 0;
    for (auto t: terms) {
        bool there_is_t = false;
        for (auto w: words)
            for (auto s: w)
                if (s.equal_symbol(t))
                    there_is_t = true;


        if (there_is_t) count_t ++;
    }*/

    /*ofstream myfile;
    myfile.open ("-conll.spice.train");
    myfile << words.size();
    myfile << " " << words.size() << endl;
    for (auto w: words) {
        if (w.size() <= 30) {
            //myfile << w.size();
            for (auto s: w)
                myfile  << s.name << " ";
            myfile << endl;
        }
    }
    myfile.close();
    exit(0);*/
    //words.erase(words.begin()+5000, words.end());
    cout << words.size() << " words" << endl;
    double avg_score = 0.0;
    double avg_log10s = 0.0;
    double avg_log2s = 0.0;
    double avg_loges = 0.0;
    double avg_log10sol = 0.0;
    double exp_avg_logs = 0.0;

    training_method = 1;
    for (int i = 0; i < 1; i ++) {
        cout << "it: " << i+1 << ": ";



        //std::random_shuffle ( words.begin(), words.end());
        vector<vector<Symbol::Symbol>> test_words;
        vector<vector<Symbol::Symbol>> train_words;
        //test_words.insert(test_words.end(), words.begin(), words.begin()+words.size()/10);
        //train_words.insert(train_words.end(), words.begin()+words.size()/10,  words.end());
        //train_words = words; //to train with all words
        if (max_word_lenght == 10) {
            train_words.insert(train_words.end(), words.begin(),  words.begin()+8751);
            test_words.insert(test_words.end(), words.begin()+8751, words.end());
        } else if (max_word_lenght == 20) {

            train_words.insert(train_words.end(), words.begin(),  words.begin()+11076);
            test_words.insert(test_words.end(), words.begin()+11076, words.end());
        }
        else if (max_word_lenght == 30){
            train_words.insert(train_words.end(), words.begin(), words.begin() +13201);
            test_words.insert(test_words.end(), words.begin() +13201, words.end());
        } else {
            train_words = words;
        }
        Grammar::Grammar g = Grammar::Grammar(terms, 3, train_words, g.pcsg, make_pair(0, 0));

        /*load_grammar(g, "music_grammar_full_beatles"+ to_string(alpha) + "_" + to_string(p_ratio)+ "_" + to_string(index_p_file) +".txt");
        //g.print_grammar();
        for (int j = 0; j < 100; j++) {
            string name = "midi_out_full_beatles";
            name += to_string(j);
            name += ".midi";
            vector<Symbol::Symbol> str = g.generate_string(30);
            for (auto &s: str)
                s = chord_terms[s.id];
            cout << "chord sequence " << j<< ": " << g.convert_vector_to_string(str) << endl << endl;
            //write_string_to_midi(str, name);
        }
        exit(0);*/
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
        cout << "start training" << endl;
        if (training_method == 0)
            g.train(g.pfa_alergia, iterations, alpha, p_ratio, time_limit);
        else if (training_method == 1)
            g.train(g.pcfg_pumping_inference, iterations, alpha, p_ratio, time_limit);
        else if (training_method == 2)
            g.train(g.pcsg_metropolis_hastings, iterations, alpha, p_ratio, time_limit);
        //g.print_grammar();
        //save_grammar(g, "cfdg_grammar_full_dyck"+ to_string(alpha) + "_" + to_string(p_ratio)+ "_" + to_string(index_p_file) +".txt");

        long double log2s = 0.0;
        long double log10s = 0.0;
        long double loges = 0.0;
        long double exp = 0.0;
        long double expI = 0.0;
        long double pcx = 0.0;
        long double log10sol = 0.0;
        for (int i2 = 0 ; i2 < test_words.size(); i2++) {
            if (training_method == 3)
                pcx = 1/ (1.0 * pow(terms.size()+1, test_words[i2].size()));
            else if (training_method == 2)
                pcx = g.find_word_probabilities_from_pcfg_inside_table(test_words[i2]);
            else
                pcx = g.find_word_probabilities(test_words[i2]);

            //cout << "word "<< i2 << " prob: "<< pcx <<" - probSol:  " << sol_pal[i2] << endl;
            //pcx = 0.0;
            if (pcx == 0.0)
                pcx = 1/ (1.0 * train_words.size() * (1.0 * pow(terms.size()+1, test_words[i2].size())));
            /*if (pcx == 0.0)
                pcx =  DBL_MIN;*/

            exp += sol_pal[i2] * log2(pcx);
            expI += sol_pal[i2] * log2(sol_pal[i2]);
            log10s -= log10(pcx);
            log2s -= log2(pcx);
            loges -= log(pcx);
            log10sol -= log10(sol_pal[i2]);

        }
        avg_score += pow(2, -exp);
        exp_avg_logs += pow(10, log10s);
        avg_log10s += log10s;
        avg_log2s += log2s;
        avg_loges += loges;
        avg_log10sol += log10sol;
        cout << " log10s: " << log10s << " - log2s: " << log2s << " - loges: " << loges<< " - logsol: " << log10sol << endl;
        cout << " s-log10s: " << pow(10,log10s) << " - s-log2s: " << pow(2,log2s) << " - s-loges: " << pow(M_E,loges)<< " - s-logsol: " << pow(10,log10sol) << endl;
        cout << " Score: " << pow(2, -exp) << " IdealScore: " << pow(2, -expI)<< endl;
        //g.print_grammar();
    }
    cout << " avg_log10s: " << avg_log10s/30 << " - avg_log2s: " << avg_log2s/30 << " - avg_loges: " << avg_loges/30<< " - avg_logsol: " << avg_log10sol/30 << endl;
    cout << " s-log10s: " << exp_avg_logs/30<< " - AVG Score: " << avg_score/30;

     //save_grammar(g, "grammar.txt");
    //Carregar Gramatica e gerar mapas
    //load_grammar(g, "cfdg_grammar_full_beatles1.500000_0.200000_21.txt");

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
    /*if (n_terminals >= 71) {
        cout << "Error. Too many terminals. Max allowed is 71 terminals" << endl;
        return;
    }*/
    for (int i = 0; i < n_terminals; i++) {
        //char symbol_name = 40;
        /*if ( index == 18)
            symbol_name+=8;*/
        if (i < 0)
            terminals.push_back(Symbol::Symbol( to_string(i), i, true, false));
        else {
            //symbol_name += i;
            //terminals.push_back(Symbol::Symbol( string(1, symbol_name), i, true, false));
            terminals.push_back(Symbol::Symbol(to_string(i), i, true, false));
        }
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
        } else if (n_symbol == 0){
            words_to_infer.push_back(word);
        }



    } while (ifs.good());
    //words_to_infer.erase(words_to_infer.end()-1, words_to_infer.end());
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
void write_string_to_midi(vector<Symbol::Symbol> string, std::string filename) {
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> starttime(0, 100);
    uniform_int_distribution<int> duration(1, 8);
    uniform_int_distribution<int> pitch(36, 84);
    uniform_int_distribution<int> velocity(40, 100);
    smf::Options options;
    options.define("n|note-count=i:10", "How many notes to randomly play");
    options.define("o|output-file=s",   "Output filename (stdout if none)");
    options.define("i|instrument=i:0",  "General MIDI instrument number");
    options.define("x|hex=b",           "Hex byte-code output");
    //options.process(argc, argv);
    smf::MidiFile midifile;
    int track   = 0;
    int channel = 0;
    int instr   = options.getInteger("instrument");
    // acounstic piano 1
    //guitar 27
    //string ensemble 49
    instr = 49;
    midifile.addTimbre(track, 0, channel, instr);
    int tpq     = midifile.getTPQ();
    int count   = options.getInteger("note-count");
    for (int i=0; i<string.size(); i++) {
        for (int j=0; j < 12; j ++) {
            if (string[i].name.compare("Other") != 0 && string[i].name[j] == '1') {
                for (int k = 0; k < 5; k++) {
                    int starttick = int((16 * (i + 1) + 4*k)/ 4.0 * tpq);
                    if (k >= 2 )
                        starttick -= 2/ 4.0 * tpq;
                    if (k >= 4 )
                        starttick -= 2/ 4.0 * tpq;
                    int key =  j+  48;

                    int endtick = starttick + int(4 / 4.0 * tpq);
                    if (k == 1 || k == 3)
                        endtick -= 2/ 4.0 * tpq;
                    midifile.addNoteOn(track, starttick, channel, key, 70);
                    midifile.addNoteOff(track, endtick, channel, key);
                }
            }
        }
    }
    midifile.sortTracks();
    if (filename.empty()) {
        if (options.getBoolean("hex")) midifile.writeHex(cout);
        else cout << midifile;
    } else
        midifile.write(filename);
}
vector<vector<Symbol::Symbol>> generate_tail_embedding_language(int max_length) {
    std::vector<Symbol::Symbol> terms = {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true,false) /*};*/, Symbol::Symbol("2", 2, true,false), Symbol::Symbol("3", 3, true,false)};
    //std::vector<Symbol::Symbol> terms = {Symbol::Symbol("V1", 0, true, false), Symbol::Symbol("V2", 1, true,false) /*};*/, Symbol::Symbol("N", 2, true,false), Symbol::Symbol("R", 3, true,false)};
    vector<vector<Symbol::Symbol>> words;
    Grammar::Grammar g = Grammar::Grammar(terms, 4, words, g.pcsg, make_pair(0, 0));
    // NT0:S NT1:NP NT2:VP NT3:CP
    // NT0 - > NT1 NT2
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right;
    std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rhs;
    rhs.first.push_back(g.non_terminals[1]);
    rhs.first.push_back(g.non_terminals[2]);
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules[0].right = right;

    // NT1 - > T2 | T2 NT3
    right.clear();
    rhs.first.clear();
    rhs.first.push_back(terms[2]);
    right.push_back(rhs);
    g.rules[1].right = right;
    rhs.first.push_back(g.non_terminals[3]);
    g.rules[1].right.push_back(rhs);

    //NT2 -> V1 | V2 NT1
    right.clear();
    rhs.first.clear();
    rhs.first.push_back(terms[0]);
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[1]);
    rhs.first.push_back(g.non_terminals[1]);
    right.push_back(rhs);
    g.rules[2].right = right;

    //NT3 -> R NT2
    right.clear();
    rhs.first.clear();
    rhs.first.push_back(terms[3]);
    rhs.first.push_back(g.non_terminals[2]);
    right.push_back(rhs);
    g.rules[3].right = right;

    g.print_grammar();
    words = g.generate_max_size_words_from_rules(max_length);
    cout << words.size() << " " << terms.size() << endl;
    for (auto w: words)
        cout << w.size() << " " << g.convert_vector_to_string(w) << endl;

    return words;

}
vector<vector<Symbol::Symbol>> generate_dyck_n_embedding_language(int n, int max_length) {
    std::vector<Symbol::Symbol> terms;
    for (int i = 0; i  < n ; i++) {
        terms.push_back(Symbol::Symbol((to_string(2*i)), 2*i, true, false));
        terms.push_back(Symbol::Symbol((to_string(2*i+1)), 2*i+1, true, false));
    }
    vector<vector<Symbol::Symbol>> words;
    Grammar::Grammar g = Grammar::Grammar(terms, 1, words, g.pcsg, make_pair(0, 0));

    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right;
    std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rhs;
    rhs.second.first = 1.0;
    g.rules[0].right = right;
    for (int i = 0; i < n; i++) {
        rhs.first.clear();
        rhs.first.push_back(terms[2*i]);
        rhs.first.push_back(g.non_terminals[0]);
        rhs.first.push_back(terms[2*i+1]);
        g.rules[0].right.push_back(rhs);

        rhs.first.clear();
        rhs.first.push_back(terms[2*i]);
        rhs.first.push_back(terms[2*i+1]);
        g.rules[0].right.push_back(rhs);
    }

    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[0]); rhs.first.push_back(g.non_terminals[0]);
    g.rules[0].right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(Symbol::Symbol("", -1, true, false));
    g.rules[0].right.push_back(rhs);
    g.print_grammar();
    words = g.generate_max_size_words_from_rules(max_length);
    cout << words.size() << " " << terms.size() << endl;
    for (auto w: words)
        cout << w.size() << " " << g.convert_vector_to_string(w) << endl;

    return words;
}
