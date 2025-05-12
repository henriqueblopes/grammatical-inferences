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
 * 25 is tail embedding (
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
 * 43 Dyck Language 6_8 with test 6_10
 * 44 Dyck Language 5_8 with test 5_10
 * 45 Dyck Language 4_8 with test 4_10
 * 46 Dyck Language 3_8 with test 3_10
 * 47 The Beatles Annotations - All Chords (144)
 * 48 The Beatles Annotations - 100 Chords
 * 49 The Beatles Annotations - 70 Chords
 * 50 The Beatles Annotations - 50 Chords
 * 51 The Beatles Annotations - 20 Chords
 * 52 Brown from NLTT 24200 - 307
 * 63-72 Brown 70-70
 * 73 |a| = |b| raw* 62103 until length 18 (should exp with 12-14)
 * 74 2|a| = |b| 3098 until 15 (should exp with 12-15)
 * 75 |a| = |b| (with 12-14) 2000 test
 * 76 2|a| = |b| (with 12-15) 2000 test*/


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
void load_spice_file(string filename, vector<Symbol::Symbol> &terminals, vector<vector<Symbol::Symbol>> &words_to_infer, int index, int n_symbol_max, bool is_numeric, vector<vector<Symbol::Symbol>> &bigger_words);
void chord_to_char(vector<Symbol::Symbol> &terminals, vector<vector<Symbol::Symbol>> &words_to_infer, vector<Symbol::Symbol> &char_symbols, vector<vector<Symbol::Symbol>> &char_words);
void write_string_to_midi(vector<Symbol::Symbol> string, std::string filename);
vector<double> load_pautomac_solution_file(string filename, int index);
vector<vector<Symbol::Symbol>> generate_palindromes(int max_length);
vector<vector<Symbol::Symbol>> generate_mod_a_eq_mod_b(int max_length);
vector<vector<Symbol::Symbol>> generate_expression_language(int max_length);
vector<vector<Symbol::Symbol>> generate_tail_embedding_language(int max_length);
vector<vector<Symbol::Symbol>> generate_dyck_n_embedding_language(int n, int max_length);
vector<vector<Symbol::Symbol>> generate_NL1(int n, int max_length);
void generate_AC1(Grammar::Grammar &g);
void generate_AC2(Grammar::Grammar &g);
void generate_NL2(Grammar::Grammar &g);
void generate_NL8(Grammar::Grammar &g);
void generate_NL8james(Grammar::Grammar &g);
/*TODO verificar P-rule Checked0 - #$% -> !#$%.189 1 1 1 1 1 1na base 68
  *
  * Experimento com spice 13 NLP (English spelling correction from Twitter Typos Corpus) (achou pumping v14 e v89)
 * And this grammar (tail embedding english)
 *  S→NP VP
 *  VP→V1 | V2 NP
 *  NP→N | N CP
 *  CP→R VP
 *
 *  */



int main(int argc, char** argv) {
    if (argc < 10)
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
    vector<vector<Symbol::Symbol>> words, chord_words, bigger_words;

    /*words = {{Symbol::Symbol("1", 1, true, false), Symbol::Symbol("1", 1, true, false), Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1, true, false)},
                                    {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("0", 0, true,false),  Symbol::Symbol("1", 1, true,false),  Symbol::Symbol("1", 1,true,false)},
                                                                            {Symbol::Symbol("0", 0, true, false), Symbol::Symbol("1", 1,true,false)}};*/
     //READ BROWN DATASET
    //NotWorking
    /*InputWords iw = InputWords(false, 20);
    iw.read_words_brown_2();
    vector<Symbol::Symbol> chord_terms = iw.generate_terminals_and_limit_string(15);
    ofstream myfile, gsMyFile, myfileMap;
    myfile.open ("52.spice.train");
    gsMyFile.open("52.spice.train.gs");
    myfileMap.open("52.spice.train.map.txt");
    myfile << chord_words.size();
    myfile << " " << chord_terms.size() << endl;
    std::random_shuffle(chord_words.begin(), chord_words.end());
    for (auto w: chord_words) {
        myfile << w.size();
        for (auto s: w) {
            myfile << " " << s.id;
            gsMyFile << s.id << " ";
        }
        myfile << endl;
        gsMyFile << endl;
    }
    for (auto s: chord_terms)
        myfileMap << s.id << " " << s.name << endl;
    exit(0);*/

    //READ MCGILL_BILLBOARDS
    /*InputWords iw = InputWords(false, 70);
    iw.read_words(false);
    //iw.input_words.erase(iw.input_words.begin(), iw.input_words.begin()+995);
    iw.iterate_chords();
    iw.change_words_to_reducted_chords();
    vector<Symbol::Symbol> chord_terms = iw.generate_terminals(iw.reducted_chord_counts);
    chord_words = iw.input_words;
    cout << chord_terms.size() << endl;
    //chord_to_char(chord_terms, chord_words, terms, words);
    ofstream myfile, gsMyFile, myfileMap;
    myfile.open ("51.spice.train");
    gsMyFile.open("51.spice.train.gs");
    myfileMap.open("51.spice.train.map.txt");
    myfile << chord_words.size();
    myfile << " " << chord_terms.size() << endl;
    //std::random_shuffle(chord_words.begin(), chord_words.end());
    for (auto w: chord_words) {
        myfile << w.size();
        for (auto s: w) {
            myfile << " " << s.id;
            gsMyFile << s.id << " ";
        }
        myfile << endl;
        gsMyFile << endl;
    }
    for (auto s: chord_terms)
        myfileMap << s.id << " " << s.name << endl;
    exit(0);*/


    //vector<Symbol::Symbol> chord_terms = iw.generate_terminals(iw.reducted_chord_counts);
     //READ MUSICAL DATABASE
    /*InputWords iw = InputWords(false, 20);
    iw.read_words_beatles(true);
    //iw.input_words.erase(iw.input_words.begin(), iw.input_words.begin()+995);
    iw.iterate_chords();
    iw.change_words_to_reducted_chords();
    vector<Symbol::Symbol> chord_terms = iw.generate_terminals(iw.reducted_chord_counts);
    chord_words = iw.input_words;
    cout << chord_terms.size() << endl;
    //chord_to_char(chord_terms, chord_words, terms, words);
    ofstream myfile, gsMyFile, myfileMap;
    myfile.open ("51.spice.train");
    gsMyFile.open("51.spice.train.gs");
    myfileMap.open("51.spice.train.map.txt");
    myfile << chord_words.size();
    myfile << " " << chord_terms.size() << endl;
    std::random_shuffle(chord_words.begin(), chord_words.end());
    for (auto w: chord_words) {
        myfile << w.size();
        for (auto s: w) {
            myfile << " " << s.id;
            gsMyFile << s.id << " ";
        }
        myfile << endl;
        gsMyFile << endl;
    }
    for (auto s: chord_terms)
        myfileMap << s.id << " " << s.name << endl;
    exit(0);*/




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
    int max_test_word_lenght = stoi((argv[9]));
    //index_p_file = 23;

    //generate_NL1(1, 11);
    //generate_NL2(1, 13);


    //SHUFFLE BROWN OR INSERT 2000 TEST
    /*load_spice_file(".spice.train", terms, words, index_p_file, max_word_lenght, false, bigger_words);
    std::random_shuffle(bigger_words.begin(), bigger_words.end());
    vector<vector<Symbol::Symbol>> bigger_limited;
    for (auto bw: bigger_words) {
        if (bw.size() <= max_test_word_lenght)
            bigger_limited.push_back(bw);
    }
    words.insert(words.end(), bigger_limited.begin(), bigger_limited.begin()+2000);
    for (int i = 1; i <= 1;i++) {
        ofstream myfile, gsMyFile, myfileMap;
        myfile.open (to_string(index_p_file)+"_2000.spice.train");
        gsMyFile.open(to_string(index_p_file)+"_2000.spice.train.gs");
        myfileMap.open(to_string(index_p_file)+"_2000.spice.train.map.txt");
        myfile << words.size();
        myfile << " " << terms.size() << endl;

        for (auto w: words) {
            myfile << w.size();
            for (auto s: w) {
                myfile << " " << s.id;
                gsMyFile << s.id << " ";
            }
            myfile << endl;
            gsMyFile << endl;
        }
        for (auto s: terms)
            myfileMap << s.id << " " << s.name << endl;
        *//*rotate(words.begin(), words.begin()+words.size()/10,words.end());*//*
    }


    exit(0);*/





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

    //CREATE DATASET BY GRAMMAR
    /*Grammar::Grammar g2 = Grammar::Grammar(terms, 1, words, g2.pcsg, std::make_pair(0, 0));
    generate_AC2(g2);
    g2.print_grammar();
    words = g2.generate_max_size_words_from_rules(max_word_lenght);
    cout << words.size() << " " << g2.terminals.size() << endl;
    for (auto w: words)
        cout << w.size() << " " << g2.convert_vector_to_string(w) << endl;
    exit(0);*/

    //load_pautomac_file(".pautomac.train", terms, words, index_p_file);
    load_spice_file(".spice.train", terms, words, index_p_file, max_word_lenght, false, bigger_words);
    vector<vector<Symbol::Symbol>> to_test;
    vector<vector<Symbol::Symbol>> to_train;
    vector<vector<Symbol::Symbol>> to_recall;
    //to_train.insert(to_train.end(), words.begin(), words.end()-1000001);
    //load_spice_file(".spice.train", terms, words, index_p_file, 19, false, bigger_words);
    to_train = words;
    for (auto bw: bigger_words) {
        if (bw.size() <= max_test_word_lenght)
            to_test.push_back(bw);
    }
    
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

    vector<double> probs;
    //CALCULAR PERP do COnn com tests especificos
    /*std::ifstream ifs ("1_pi_6_8_probs.txt", std::ifstream::in);
    std::ifstream ifs2 ("2_pi_6_8_probs.txt", std::ifstream::in);*/
    /*std::ifstream ifs ("experimento208.txt", std::ifstream::in);
    //std::ifstream ifs2 ("jpcfg-probs2-20.txt", std::ifstream::in);
    std::ifstream ifs2 ("1_3_rand_tr_conll_2_20.txt", std::ifstream::in);
    string line;
    string line2;
    double loglike = 0.0;
    double loglike2 = 0.0;
    double negloglike = 0.0;
    double negloglike2 = 0.0;
    double n_words = 0.0;
    double n_words_pos = 0.0;
    double n_words_neg = 0.0;
    while (ifs.good()) {
        n_words+=1.0;
        getline(ifs,line);
        getline(ifs2,line2);
        if (line.find(" 0 -") == string::npos) {
            n_words_pos += 1.0;
            string aux = line.substr(line.find(":")+1, line.find(" - probSol:")-line.find(":")-1);
            string aux2 = line2.substr(line2.find(":")+1, line2.find(" - probSol:")-line2.find(":")-1);
            probs.push_back(stod(aux));
            if (!aux.empty()) {
                cout << "logPI: " << -log10(stod(aux));
                //cout << " - logPcfg: " <<  -log10(stod(aux2)) << endl;
                cout << endl;
                loglike -= log10(stod(aux));
                //loglike2 -= log10(stod(aux2));

            }
        } else {
            n_words_neg += 1.0;
            string aux = line.substr(line.find(":")+1, line.find(" - probSol:")-line.find(":")-1);
            string aux2 = line2.substr(line2.find(":")+1, line2.find(" - probSol:")-line2.find(":")-1);
            probs.push_back(stod(aux));
            if ( !aux2.empty()) {
                cout << "logPI: " << -log10(stod(aux));
                //cout << " - logPcfg: " <<  -log10(stod(aux2)) << endl;
                cout << endl;
                negloglike -= log10(stod(aux));
                //negloglike2 -= log10(stod(aux2));
            }
        }
        if (n_words >=1000000.0)
            break;
    }*/
    //cout << "logPi: " << loglike+negloglike2<< " logPiPos: " << loglike//*(n_words_pos-1) <<" - logpcfg: " << loglike2//*(2428-n_words_neg) <<  " - negloglike: " << negloglike//*(2428-n_words_neg) << " - negloglike2: " << negloglike2//*(2428-n_words_neg) << " n_words: " << n_words-1 << " n_words_pos: " << n_words_pos-1 << " n_words_neg: " << n_words_neg<<endl;
    /*return 0;*/
    /*for (int i = 51; i <= 53; i++) {
        terms.clear();
        words.clear();
        load_spice_file(".spice.train", terms, words, i, max_word_lenght);
        std::random_shuffle ( words.begin(), words.end());
        ofstream myfile, myfile2;
        myfile.open (to_string(i)+".dyck.spice.train");
        myfile << words.size();
        myfile << " " << terms.size() << endl;
        for (auto w: words) {
            if (w.size() <= 30) {
                myfile << w.size() ;
                for (auto s: w)
                    myfile  << " " << s.name;
                myfile << endl;
            }
        }
        myfile.close();

        myfile2.open (to_string(i-50)+"-dyck.pcfg-dyck-5-10.yld");
        for (auto w: words) {
            if (w.size() <= 30) {
                for (auto s: w)
                    myfile2  << s.name << " ";
                myfile2 << endl;
            }
        }
        myfile2.close();
    }
    exit(0);*/
    //words.erase(words.begin()+5000, words.end());
    cout << words.size() << " words" << endl;
    double avg_score = 0.0;
    double avg_log10s = 0.0;
    double avg_log2s = 0.0;
    double avg_loges = 0.0;
    double avg_log10sol = 0.0;
    double exp_avg_logs = 0.0;

    //training_method = 1;
    for (int i = 0; i < 1; i ++) {
        cout << "it: " << i+1 << ": ";



        //std::random_shuffle ( words.begin(), words.end());


        vector<vector<Symbol::Symbol>> test_words;
        vector<vector<Symbol::Symbol>> train_words;
        //Dyck_6_8
        /*train_words.insert(train_words.end(), words.begin(),  words.begin()+19303);
        test_words.insert(test_words.end(), words.begin()+19303, words.end());*/

        //Dyck_5_8
        /*train_words.insert(train_words.end(), words.begin(),  words.begin()+9431);
        test_words.insert(test_words.end(), words.begin()+9431, words.end());*/

        //Dyck_4_8
        /*train_words.insert(train_words.end(), words.begin(),  words.begin()+3941);
        test_words.insert(test_words.end(), words.begin()+3941, words.end());*/

        //Dyck_3_8
        /*train_words.insert(train_words.end(), words.begin(),  words.begin()+1291);
        test_words.insert(test_words.end(), words.begin()+1291, words.end());*/

        /*test_words.insert(test_words.end(), words.begin(), words.begin()+330913);
        train_words.insert(train_words.end(), words.begin()+330913,  words.end());*/

        /*test_words.insert(test_words.end(), words.begin(), words.begin()+words.size()/10);
        train_words.insert(train_words.end(), words.begin()+words.size()/10,  words.end());*/


        train_words = words; //to train with all words


        //train_words.erase(train_words.begin()+1000000, train_words.end());
        //test_words = bigger_words;

        test_words = to_test;
        train_words = to_train;
        /*if (max_word_lenght == 10) {
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
        }*/
        train_words.push_back(vector<Symbol::Symbol>());
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
        else if (training_method == 4)
            g.g_tp = g.pcfg;
        else
            exit(-1);

        unordered_map<string, int> sol_count;
        unordered_map<string, double> sol_pal;
        for (auto w: test_words) {
            if (sol_count.find(g.convert_vector_to_string(w)) == 0)
                sol_count[g.convert_vector_to_string(w)] = 1;
            else
                sol_count[g.convert_vector_to_string(w)] ++;
        }
        /*for (auto w: test_words) {
            sol_pal[g.convert_vector_to_string(w)] = sol_count[g.convert_vector_to_string(w)]*1.0/1000000.0;
        }*/


        cout << "start training" << endl;
        if (training_method == 0)
            g.train(g.pfa_alergia, iterations, alpha, p_ratio, time_limit);
        else if (training_method == 1)
            g.train(g.pcfg_pumping_inference, iterations, alpha, p_ratio, time_limit);
        else if (training_method == 2)
            g.train(g.pcsg_metropolis_hastings, iterations, alpha, p_ratio, time_limit);
        else if (training_method == 4)
            //load_grammar(g, "pcfg_grammar_" + to_string(index_p_file) + "_" + to_string(alpha) + "_" + to_string(p_ratio)+ "_" + to_string(index_p_file) +"_DEL.txt");
            cout << "training complete" << endl;
        g.print_grammar();
        if (training_method == 1)
            save_grammar(g, "pcfg_grammar_"+ to_string(index_p_file) + "_" + to_string(alpha) + "_" + to_string(p_ratio)+ "_" + to_string(index_p_file) +".txt");
        //exit(1);



        double recall = 0.0;
        Grammar::Grammar g2 = Grammar::Grammar(terms, 1, words, g.pcsg, std::make_pair(0, 0));
        generate_AC1(g2);
        g2.n_non_terminals = g2.non_terminals.size();
        g2.print_grammar();
        g2.convert_to_cnf_full();
        g2.print_grammar();
        for (int j = 0 ; j < to_test.size(); j++) {
            vector<Symbol::Symbol> word_recall = g.generate_string(20);
            cout << "Word Recall "<< j <<": " << g2.convert_vector_to_string(word_recall);
            double wr_prob = g2.probabilistic_cky(word_recall);
            cout <<" - Prob: " << wr_prob << endl;
            if (wr_prob >0.0)
                recall +=1.0;
        }
        recall /= to_test.size();

        /*Grammar::Grammar g2 = Grammar::Grammar(terms, 3, train_words, g.pcsg, make_pair(0, 0));;
        load_grammar(g2, "pcfg_grammar_" + to_string(index_p_file) + "_" + to_string(alpha) + "_" + to_string(p_ratio)+ "_" + to_string(index_p_file) +".txt");
        g2.rules = g.rules;
        g2.n_non_terminals = g.non_terminals.size();
        g2.n_terminals = g.terminals.size();
        g2.convert_to_cnf_full();*/
        g.convert_to_cnf_full();
        /*Grammar::Grammar g2 = Grammar::Grammar(terms, 1, words, g.pcsg, std::make_pair(0, 0));
        generate_NL8(g2);
        Grammar::Grammar g3 = Grammar::Grammar(terms, 1, words, g.pcsg, std::make_pair(0, 0));
        generate_NL8james(g3);
        g2.n_non_terminals = g2.non_terminals.size();
        g3.n_non_terminals = g3.non_terminals.size();
        g2.convert_to_cnf_full();
        g3.convert_to_cnf_full();
        g2.print_grammar();
        g3.print_grammar();*/
        /*if (training_method == 1)
            save_grammar(g, "pcfg_grammar_"+ to_string(index_p_file) + "_" + to_string(alpha) + "_" + to_string(p_ratio)+ "_" + to_string(index_p_file) +"_DEL.txt");
        exit(1);*/



        long double log2s_add1_smooth = 0.0;
        long double log10s_add1_smooth = 0.0;
        long double loges_add1_smooth = 0.0;
        long double log2s_random = 0.0;
        long double log10s_random = 0.0;
        long double loges_random = 0.0;
        long double log2s_worse_random = 0.0;
        long double log10s_worse_random = 0.0;
        long double loges_worse_random = 0.0;
        long double log2s_randsol = 0.0;
        long double log10s_randsol = 0.0;
        long double loges_randsol = 0.0;
        long double log2s_Wrandsol = 0.0;
        long double log10s_Wrandsol = 0.0;
        long double loges_Wrandsol = 0.0;
        long double log2s_Add1sol = 0.0;
        long double log10s_Add1sol = 0.0;
        long double loges_Add1sol = 0.0;
        long double exp = 0.0;
        long double expI = 0.0;
        long double pcx = 0.0;
        long double logesol = 0.0;
        long double log2sol = 0.0;
        long double log10sol = 0.0;
        long double kullbackD = 0.0;
        long double kullback_add1_smooth = 0.0;
        long double kullback_random = 0.0;
        long double kullback_worse_random = 0.0;
        long double kullback_randsol = 0.0;
        long double kullback_Wrandsol = 0.0;
        long double kullback_Add1sol = 0.0;
        long double amount_zeros = 0.0;
        long double amount_bigger_length = 0.0;
        long double amount_near_01_k = 0.0;
        for (int i2 = 0 ; i2 < test_words.size(); i2++) {
            if  (test_words[i2].size() > max_test_word_lenght) {
                cout << "word " << i2 << ": length bigger than " << max_test_word_lenght << endl;
                amount_bigger_length+=1.0;
            } else {
                if (training_method == 4) {
                    pcx = 0.0;
                    //pcx = probs[i2];
                    //pcx = g3.probabilistic_cky(test_words[i2]);
                }
                else if (training_method == 3)
                    pcx = 1/ (1.0 * pow(terms.size()+1, test_words[i2].size()));
                else if (training_method == 2)
                    pcx = g.find_word_probabilities_from_pcfg_inside_table(test_words[i2]);
                else {
                    //pcx = g.find_word_probabilities(test_words[i2]);
                    //pcx = g2.probabilistic_cky(test_words[i2]);
                    pcx = g.probabilistic_cky(test_words[i2]);
                    //probs[i2] = pcx;
                }
                double pcx2 = 1.0/(1.0*test_words.size());//

                //pcx2 =  g2.probabilistic_cky(test_words[i2]);
                //cout << pcx2 << endl;
                sol_pal[g.convert_vector_to_string(test_words[i2])] = pcx2;
                double add1 = (test_words.size()+ sol_count[g.convert_vector_to_string(test_words[i2])]);
                double pcx_add1_smooth, pcx_random, pcx_worse_random;
                double random_g = 1/ (1.0 *(1.0 * pow(terms.size()+1, test_words[i2].size())));
                double wrandom_g = 1/ (1.0 * train_words.size() * (1.0 * pow(terms.size()+1, test_words[i2].size())));
                //cout << "word "<< i2+1 << " prob: "<< probs[i2] << " - probSol:  " << random_g /*<< " probcky: " << g2.probabilistic_cky(test_words[i2])*/ << endl;

                //pcx = 0.0;

                //kullbackD += log(random_g/add1) * random_g;
                pcx_add1_smooth = (pcx + 1)/(1.0*(test_words.size()+ sol_count[g.convert_vector_to_string(test_words[i2])]));
                pcx_worse_random = 1/ (1.0 * train_words.size() * (1.0 * pow(terms.size()+1, test_words[i2].size())));
                pcx_random =       1/ (1.0 *(1.0 * pow(terms.size()+1, test_words[i2].size())));

                cout << "word "<< i2+1 << " - prob: "<<  pcx << " - probckyadd1: " << pcx_add1_smooth << " - probSolRand:  " << pcx_random << " - probSolWRand: " << pcx_worse_random << " - probSolTest:  " << sol_pal[g.convert_vector_to_string(test_words[i2])]  << endl;
                if (pcx == 0.0) {
                    //PCX igual a add1 ou aleatorio ou pior que aleatorio
                    amount_zeros += 1.0;

                    pcx = pcx_add1_smooth;
                    kullback_add1_smooth +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                    log10s_add1_smooth -= log10(pcx_add1_smooth);
                    log2s_add1_smooth -= log2(pcx_add1_smooth);
                    loges_add1_smooth -= log(pcx_add1_smooth);

                    pcx = pcx_random;
                    kullback_random +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                    log10s_random -= log10(pcx_random);
                    log2s_random -= log2(pcx_random);
                    loges_random -= log(pcx_random);

                    pcx = pcx_worse_random;
                    kullback_worse_random +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                    log10s_worse_random -= log10(pcx_worse_random);
                    log2s_worse_random -= log2(pcx_worse_random);
                    loges_worse_random -= log(pcx_worse_random);


                } else {
                    kullback_add1_smooth +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                    kullback_random +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                    kullback_worse_random +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                    log10s_add1_smooth -= log10(pcx);
                    log2s_add1_smooth -= log2(pcx);
                    loges_add1_smooth -= log(pcx);
                    log10s_random -= log10(pcx);
                    log2s_random -= log2(pcx);
                    loges_random -= log(pcx);
                    log10s_worse_random -= log10(pcx);
                    log2s_worse_random -= log2(pcx);
                    loges_worse_random -= log(pcx);
                }

                log10sol -= log10(sol_pal[g.convert_vector_to_string(test_words[i2])]);
                log2sol -= log2(sol_pal[g.convert_vector_to_string(test_words[i2])]);
                logesol -= log(sol_pal[g.convert_vector_to_string(test_words[i2])]);

                kullback_Add1sol +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx_add1_smooth)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                log10s_Add1sol -= log10(pcx_add1_smooth);
                log2s_Add1sol -= log2(pcx_add1_smooth);
                loges_Add1sol -= log(pcx_add1_smooth);

                kullback_randsol +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx_random)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                log10s_randsol -= log10(pcx_random);
                log2s_randsol -= log2(pcx_random);
                loges_randsol -= log(pcx_random);

                kullback_Wrandsol +=  log(sol_pal[g.convert_vector_to_string(test_words[i2])]/pcx_worse_random)  * sol_pal[g.convert_vector_to_string(test_words[i2])];
                log10s_Wrandsol -= log10(pcx_worse_random);
                log2s_Wrandsol -= log2(pcx_worse_random);
                loges_Wrandsol -= log(pcx_worse_random);
            }
        }



        cout << "My Sol: " << endl;
        cout << " KB_add1: " << kullback_add1_smooth <<  " - KB_random: " << kullback_random << " - KB_Wrandom: "<< kullback_worse_random << endl;
        cout << " log10sAVGAdd1: " << log10s_add1_smooth /test_words.size() << " - log2sAVGAdd1: " << log2s_add1_smooth /test_words.size() << " - logesAVGAdd1: " << loges_add1_smooth /test_words.size()<< endl;
        cout << " log10s_AVGrandom: " << log10s_random/test_words.size() << " - log2s_AVGrandom: " << log2s_random/test_words.size() << " - loges_AVGrandom: " << loges_random/test_words.size() << endl;
        cout << " log10s_AVGWrandom: " << log10s_worse_random/test_words.size() << " - log2s_AVGWrandom: " << log2s_worse_random/test_words.size() << " - loges_AVGWrandom: " << loges_worse_random/test_words.size() << endl << endl;

        cout << " Test data size: " << test_words.size() << endl;
        cout << " KB_add1: " << kullback_add1_smooth <<  " - KB_random: " << kullback_random << " - KB_Wrandom: "<< kullback_worse_random << endl;
        cout << " KB_add1sol: " << kullback_Add1sol <<  " - KB_randsol: " << kullback_randsol << " - KB_Wrandsol: "<< kullback_Wrandsol << endl << endl;

        cout << " log10sAVGAdd1: " << log10s_add1_smooth /test_words.size() << " - log2sAVGAdd1: " << log2s_add1_smooth /test_words.size() << " - logesAVGAdd1: " << loges_add1_smooth /test_words.size()<< endl;
        cout << " log10s_AVGrandom: " << log10s_random/test_words.size() << " - log2s_AVGrandom: " << log2s_random/test_words.size() << " - loges_AVGrandom: " << loges_random/test_words.size() << endl;
        cout << " log10s_AVGWrandom: " << log10s_worse_random/test_words.size() << " - log2s_AVGWrandom: " << log2s_worse_random/test_words.size() << " - loges_AVGWrandom: " << loges_worse_random/test_words.size() << endl << endl;
        cout << " log10sAVGAdd1SOL: " << log10s_Add1sol /test_words.size() << " - log2sAVGAdd1SOL: " << log2s_Add1sol /test_words.size() << " - logesAVGAdd1SOL: " << loges_Add1sol /test_words.size()<<  endl;
        cout << " log10s_AVGrandomSOL: " << log10s_randsol/test_words.size() << " - log2s_AVGrandomSOL: " << log2s_randsol/test_words.size() << " - loges_AVGrandomSOL: " << loges_randsol/test_words.size() << endl;
        cout << " log10sol: " << log10sol/test_words.size() << " - log2sol: " << log2sol/test_words.size() << " - logesol: " << logesol/test_words.size()  << endl;
        cout << " Zero Probs: " << amount_zeros << " - Precision: " << (test_words.size()-amount_bigger_length-amount_zeros)/(test_words.size()-amount_bigger_length) << " - Recall: " << recall << endl;
        //g.print_grammar();
    }

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
    g.n_non_terminals = g.non_terminals.size();
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

void load_spice_file(string filename, vector<Symbol::Symbol> &terminals, vector<vector<Symbol::Symbol>> &words_to_infer, int index, int n_symbol_max, bool is_numeric, vector<vector<Symbol::Symbol>> &bigger_words) {
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
    if (is_numeric){
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
    }
    do  {
        getline(ifs,line);
        size_t tokenPos = line.find(" ");
        int n_symbol = 0;
        if (tokenPos != string::npos)
            n_symbol = stoi(line.substr(0,tokenPos+1));
        vector<Symbol::Symbol> word;

        if ( n_symbol > 0 ) {
            line = line.substr(tokenPos+1, string::npos);
            while (tokenPos != string::npos) {

                tokenPos = line.substr(0, string::npos).find(" ");
                int symbol = 0;
                if (is_numeric)
                    int symbol = stoi(line.substr(0,tokenPos));
                else {
                    bool exist_s = false;
                    for (auto s2: terminals) {
                        if (s2.name.compare(line.substr(0,tokenPos)) == 0) {
                            exist_s = true;
                            symbol = s2.id;
                            break;
                        }
                    }
                    if (!exist_s) {
                        Symbol::Symbol s = Symbol::Symbol(line.substr(0,tokenPos), terminals.size(), true, false);
                        terminals.push_back(s);
                        symbol = s.id;
                    }

                }
                word.push_back(terminals[symbol]);
                line = line.substr(tokenPos+1, string::npos);
            }
            if (word.size() <= n_symbol_max)
                words_to_infer.push_back(word);
            else
                bigger_words.push_back(word);
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

vector<vector<Symbol::Symbol>> generate_NL1(int n, int max_length) {
    std::vector<Symbol::Symbol> terms;
    terms.push_back(Symbol::Symbol("the", 0, true, false));
    terms.push_back(Symbol::Symbol("ate", 1, true, false));
    terms.push_back(Symbol::Symbol("slept", 2, true, false));
    terms.push_back(Symbol::Symbol("saw", 3, true, false));
    terms.push_back(Symbol::Symbol("heard", 4, true, false));
    terms.push_back(Symbol::Symbol("cat", 5, true, false));
    terms.push_back(Symbol::Symbol("dog", 6, true, false));
    terms.push_back(Symbol::Symbol("big", 7, true, false));
    terms.push_back(Symbol::Symbol("old", 8, true, false));

    vector<vector<Symbol::Symbol>> words;
    Grammar::Grammar g = Grammar::Grammar(terms, 1, words, g.pcsg, make_pair(0, 0));
    g.non_terminals.clear();
    g.rules.clear();

    g.non_terminals.push_back(Symbol::Symbol("S", 0, false, false));
    g.non_terminals.push_back(Symbol::Symbol("VP", 1, false, false));
    g.non_terminals.push_back(Symbol::Symbol("NP", 2, false, false));
    g.non_terminals.push_back(Symbol::Symbol("AP", 3, false, false));
    g.non_terminals.push_back(Symbol::Symbol("VERBI", 4, false, false));
    g.non_terminals.push_back(Symbol::Symbol("VERBT", 5, false, false));
    g.non_terminals.push_back(Symbol::Symbol("NOUN", 6, false, false));
    g.non_terminals.push_back(Symbol::Symbol("ADJ", 7, false, false));

    vector<Symbol::Symbol> left;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;

    // S -> NP VP
    left.push_back(g.non_terminals[0]);  // S
    rhs.first.push_back(g.non_terminals[2]);  // NP
    rhs.first.push_back(g.non_terminals[1]);  // VP
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // VP -> VERBI
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[1]);  // VP
    rhs.first.push_back(g.non_terminals[4]);  // VERBI
    right.push_back(rhs);
    // VP -> VERBT NP
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[5]);  // VERBT
    rhs.first.push_back(g.non_terminals[2]);  // NP
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // NP -> 'the' NOUN
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[2]);  // NP
    rhs.first.push_back(terms[0]);  // 'the'
    rhs.first.push_back(g.non_terminals[6]);  // NOUN
    right.push_back(rhs);

    // NP -> 'the' AP NOUN
    rhs.first.clear();
    rhs.first.push_back(terms[0]);  // 'the'
    rhs.first.push_back(g.non_terminals[3]);  // AP
    rhs.first.push_back(g.non_terminals[6]);  // NOUN
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // AP -> ADJ
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[3]);  // AP
    rhs.first.push_back(g.non_terminals[7]);  // ADJ
    right.push_back(rhs);
    // AP -> ADJ AP
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[7]);  // ADJ
    rhs.first.push_back(g.non_terminals[3]);  // AP
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // VERBI -> 'ate' | 'slept'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[4]);  // VERBI
    rhs.first.push_back(terms[1]);  // 'ate'
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[2]);  // 'slept'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // VERBT -> 'saw' | 'heard'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[5]);  // VERBT
    rhs.first.push_back(terms[3]);  // 'saw'
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[4]);  // 'heard'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // NOUN -> 'cat' | 'dog'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[6]);  // NOUN
    rhs.first.push_back(terms[5]);  // 'cat'
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[6]);  // 'dog'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // ADJ -> 'big' | 'old'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[7]);  // ADJ
    rhs.first.push_back(terms[7]);  // 'big'
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[8]);  // 'old'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    g.print_grammar();

    std::ofstream outFile("/home/petcomp/CLionProjects/grammatical-inferences/examples/palavrasNL1_" + to_string(max_length) + ".txt");
    words = g.generate_max_size_words_from_rules(max_length);
    cout << words.size() << " " << terms.size() << endl;
    outFile << words.size() << " " << terms.size() << std::endl;
    for (auto w: words) {
        cout << w.size() << " " << g.convert_vector_to_string(w) << endl;
        outFile << w.size() << " " << g.convert_vector_to_string(w) << std::endl;
    }
    outFile.close();

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
    //g.print_grammar();
    words = g.generate_max_size_words_from_rules(max_length);
    cout << words.size() << " " << terms.size() << endl;
    for (auto w: words)
        cout << w.size() << " " << g.convert_vector_to_string(w) << endl;

    return words;
}

void generate_NL8(Grammar::Grammar &g) {
    std::vector<Symbol::Symbol> terms;

    // Terminais
    terms.push_back(Symbol::Symbol("the", 0, true, false));
    terms.push_back(Symbol::Symbol("car", 1, true, false));
    terms.push_back(Symbol::Symbol("city", 2, true, false));
    terms.push_back(Symbol::Symbol("house", 3, true, false));
    terms.push_back(Symbol::Symbol("shop", 4, true, false));
    terms.push_back(Symbol::Symbol("large", 5, true, false));
    terms.push_back(Symbol::Symbol("small", 6, true, false));
    terms.push_back(Symbol::Symbol("ugly", 7, true, false));
    terms.push_back(Symbol::Symbol("beautiful", 8, true, false));
    terms.push_back(Symbol::Symbol("knows", 9, true, false));
    terms.push_back(Symbol::Symbol("likes", 10, true, false));
    terms.push_back(Symbol::Symbol("misses", 11, true, false));
    terms.push_back(Symbol::Symbol("sees", 12, true, false));
    terms.push_back(Symbol::Symbol("thinks", 13, true, false));
    terms.push_back(Symbol::Symbol("hopes", 14, true, false));
    terms.push_back(Symbol::Symbol("tells", 15, true, false));
    terms.push_back(Symbol::Symbol("says", 16, true, false));
    terms.push_back(Symbol::Symbol("appears", 17, true, false));
    terms.push_back(Symbol::Symbol("is", 18, true, false));
    terms.push_back(Symbol::Symbol("seems", 19, true, false));
    terms.push_back(Symbol::Symbol("looks", 20, true, false));
    terms.push_back(Symbol::Symbol("with", 21, true, false));
    terms.push_back(Symbol::Symbol("near", 22, true, false));
    terms.push_back(Symbol::Symbol("in", 23, true, false));
    terms.push_back(Symbol::Symbol("from", 24, true, false));
    terms.push_back(Symbol::Symbol("John", 25, true, false));
    terms.push_back(Symbol::Symbol("Mary", 26, true, false));
    terms.push_back(Symbol::Symbol("the", 27, true, false));
    terms.push_back(Symbol::Symbol("man", 28, true, false));
    terms.push_back(Symbol::Symbol("child", 29, true, false));
    terms.push_back(Symbol::Symbol("that", 30, true, false));

    std::vector<std::vector<Symbol::Symbol>> words;
    g.non_terminals.clear();
    g.rules.clear();

    // Não-terminais
    g.non_terminals.push_back(Symbol::Symbol("S", 0, false, false));
    g.non_terminals.push_back(Symbol::Symbol("NP", 1, false, false));
    g.non_terminals.push_back(Symbol::Symbol("Vi", 2, false, false));
    g.non_terminals.push_back(Symbol::Symbol("ADV", 3, false, false));
    g.non_terminals.push_back(Symbol::Symbol("NPa", 4, false, false));
    g.non_terminals.push_back(Symbol::Symbol("VPa", 5, false, false));
    g.non_terminals.push_back(Symbol::Symbol("Vs", 6, false, false));
    g.non_terminals.push_back(Symbol::Symbol("NPp", 7, false, false));
    g.non_terminals.push_back(Symbol::Symbol("Vt", 8, false, false));
    g.non_terminals.push_back(Symbol::Symbol("P", 9, false, false));

    std::vector<Symbol::Symbol> left;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;

    // S -> NP Vi ADV
    left.push_back(g.non_terminals[0]); // S
    rhs.first.push_back(g.non_terminals[1]); // NP
    rhs.first.push_back(g.non_terminals[2]); // Vi
    rhs.first.push_back(g.non_terminals[3]); // ADV
    rhs.second.first = 1.0;
    right.push_back(rhs);

    // S -> NPa VPa
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[4]); // NPa
    rhs.first.push_back(g.non_terminals[5]); // VPa
    right.push_back(rhs);

    // S -> NPa Vs 'that' S
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[4]); // NPa
    rhs.first.push_back(g.non_terminals[6]); // Vs
    rhs.first.push_back(terms[30]);          // 'that'
    rhs.first.push_back(g.non_terminals[0]); // S
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // NP -> NPa
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[1]); // NP
    rhs.first.push_back(g.non_terminals[4]); // NPa
    right.push_back(rhs);

    // NP -> NPp
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[7]); // NPp
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // Vi -> 'appears' | 'is' | 'seems' | 'looks'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[2]); // Vi
    rhs.first.push_back(terms[17]); // 'appears'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[18]); // 'is'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[19]); // 'seems'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[20]); // 'looks'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // ADV -> 'large' | 'small' | 'ugly' | 'beautiful'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[3]); // ADV
    rhs.first.push_back(terms[5]); // 'large'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[6]); // 'small'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[7]); // 'ugly'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[8]); // 'beautiful'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // NPa -> 'John' | 'Mary' | 'the' 'man' | 'the' 'child'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[4]); // NPa
    rhs.first.push_back(terms[25]); // 'John'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[26]); // 'Mary'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[27]); // 'the'
    rhs.first.push_back(terms[28]); // 'man'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[27]); // 'the'
    rhs.first.push_back(terms[29]); // 'child'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // VPa -> Vt NP
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[5]); // VPa
    rhs.first.push_back(g.non_terminals[8]); // Vt
    rhs.first.push_back(g.non_terminals[1]); // NP
    right.push_back(rhs);

    // VPa -> Vt NP P NPp
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[8]); // Vt
    rhs.first.push_back(g.non_terminals[1]); // NP
    rhs.first.push_back(g.non_terminals[9]); // P
    rhs.first.push_back(g.non_terminals[7]); // NPp
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // Vs -> 'thinks' | 'hopes' | 'tells' | 'says'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[6]); // Vs
    rhs.first.push_back(terms[13]); // 'thinks'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[14]); // 'hopes'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[15]); // 'tells'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[16]); // 'says'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));


    // NPp -> 'the' 'car' | 'the' 'city' | 'the' 'house' | 'the' 'shop'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[7]); // NPp
    rhs.first.push_back(terms[0]); // 'the'
    rhs.first.push_back(terms[1]); // 'car'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[0]); // 'the'
    rhs.first.push_back(terms[2]); // 'city'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[0]); // 'the'
    rhs.first.push_back(terms[3]); // 'house'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[0]); // 'the'
    rhs.first.push_back(terms[4]); // 'shop'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));



    // Vt -> 'knows' | 'likes' | 'misses' | 'sees'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[8]); // Vt
    rhs.first.push_back(terms[9]); // 'knows'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[10]); // 'likes'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[11]); // 'misses'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[12]); // 'sees'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));




    // P -> 'with' | 'near' | 'in' | 'from'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[9]); // P
    rhs.first.push_back(terms[21]); // 'with'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[22]); // 'near'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[23]); // 'in'
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[24]); // 'from'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    g.normalize_probs();
    // Imprime a gramática
    g.terminals = terms;
    g.print_grammar();
    /*std::ofstream outFile("/home/petcomp/CLionProjects/grammatical-inferences/examples/palavrasNL8_" + to_string(max_length) + ".txt");
    words = g.generate_max_size_words_from_rules(max_length);
    cout << words.size() << " " << terms.size() << endl;
    outFile << words.size() << " " << terms.size() << std::endl;
    for (auto w: words) {
        cout << w.size() << " " << g.convert_vector_to_string(w) << endl;
        outFile << w.size() << " " << g.convert_vector_to_string(w) << std::endl;
    }
    outFile.close();
    return words;*/
}

void generate_NL2(Grammar::Grammar &g) {
    std::vector<Symbol::Symbol> terms;
    terms.push_back(Symbol::Symbol("the", 0, true, false));
    terms.push_back(Symbol::Symbol("a", 1, true, false));
    terms.push_back(Symbol::Symbol("heard", 2, true, false));
    terms.push_back(Symbol::Symbol("saw", 3, true, false));
    terms.push_back(Symbol::Symbol("cat", 4, true, false));
    terms.push_back(Symbol::Symbol("dog", 5, true, false));
    terms.push_back(Symbol::Symbol("mouse", 6, true, false));
    terms.push_back(Symbol::Symbol("that", 7, true, false));

    vector<vector<Symbol::Symbol>> words;
    g.non_terminals.clear();
    g.rules.clear();
    g.non_terminals.push_back(Symbol::Symbol("S", 0, false, false));
    g.non_terminals.push_back(Symbol::Symbol("VP", 1, false, false));
    g.non_terminals.push_back(Symbol::Symbol("NP", 2, false, false));
    g.non_terminals.push_back(Symbol::Symbol("RC", 3, false, false));
    g.non_terminals.push_back(Symbol::Symbol("VERB", 4, false, false));
    g.non_terminals.push_back(Symbol::Symbol("NOUN", 5, false, false));
    g.non_terminals.push_back(Symbol::Symbol("ART", 6, false, false));
    g.non_terminals.push_back(Symbol::Symbol("REL", 7, false, false));

    vector<Symbol::Symbol> left;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;

    // S -> NP VP
    left.push_back(g.non_terminals[0]);  // S
    rhs.first.push_back(g.non_terminals[2]);  // NP
    rhs.first.push_back(g.non_terminals[1]);  // VP
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // VP -> VERB NP
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[1]);  // VP
    rhs.first.push_back(g.non_terminals[4]);  // VERB
    rhs.first.push_back(g.non_terminals[2]);  // NP
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // NP -> ART NOUN
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[2]);  // NP
    rhs.first.push_back(g.non_terminals[6]);  // ART
    rhs.first.push_back(g.non_terminals[5]);  // NOUN
    right.push_back(rhs);
    // NP -> ART NOUN RC
    rhs.first.push_back(g.non_terminals[3]);  // RC
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // RC -> REL VP
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[3]);  // RC
    rhs.first.push_back(g.non_terminals[7]);  // REL
    rhs.first.push_back(g.non_terminals[1]);  // VP
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // VERB -> 'heard' | 'saw'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[4]);  // VERB
    rhs.first.push_back(terms[2]);  // 'heard'
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[3]);
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // NOUN -> 'cat' | 'dog' | 'mouse'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[5]);  // NOUN
    rhs.first.push_back(terms[4]);  // 'cat'
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[5]);  // 'dog'
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[6]);  // 'mouse'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // ART -> 'the' | 'a'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[6]);  // ART
    rhs.first.push_back(terms[0]);  // 'the'
    right.push_back(rhs);
    rhs.first.clear();
    rhs.first.push_back(terms[1]);  // 'a'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // REL -> 'that'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[7]);  // REL
    rhs.first.push_back(terms[7]);  // 'that'
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // Imprime a gramática
    g.print_grammar();


    g.normalize_probs();
    // Imprime a gramática
    g.terminals = terms;
    g.print_grammar();
}


void generate_NL8james(Grammar::Grammar &g) {
    std::vector<Symbol::Symbol> terms;

    // Terminais
    terms.push_back(Symbol::Symbol("the", 0, true, false));
    terms.push_back(Symbol::Symbol("car", 1, true, false));
    terms.push_back(Symbol::Symbol("city", 2, true, false));
    terms.push_back(Symbol::Symbol("house", 3, true, false));
    terms.push_back(Symbol::Symbol("shop", 4, true, false));
    terms.push_back(Symbol::Symbol("large", 5, true, false));
    terms.push_back(Symbol::Symbol("small", 6, true, false));
    terms.push_back(Symbol::Symbol("ugly", 7, true, false));
    terms.push_back(Symbol::Symbol("beautiful", 8, true, false));
    terms.push_back(Symbol::Symbol("knows", 9, true, false));
    terms.push_back(Symbol::Symbol("likes", 10, true, false));
    terms.push_back(Symbol::Symbol("misses", 11, true, false));
    terms.push_back(Symbol::Symbol("sees", 12, true, false));
    terms.push_back(Symbol::Symbol("thinks", 13, true, false));
    terms.push_back(Symbol::Symbol("hopes", 14, true, false));
    terms.push_back(Symbol::Symbol("tells", 15, true, false));
    terms.push_back(Symbol::Symbol("says", 16, true, false));
    terms.push_back(Symbol::Symbol("appears", 17, true, false));
    terms.push_back(Symbol::Symbol("is", 18, true, false));
    terms.push_back(Symbol::Symbol("seems", 19, true, false));
    terms.push_back(Symbol::Symbol("looks", 20, true, false));
    terms.push_back(Symbol::Symbol("with", 21, true, false));
    terms.push_back(Symbol::Symbol("near", 22, true, false));
    terms.push_back(Symbol::Symbol("in", 23, true, false));
    terms.push_back(Symbol::Symbol("from", 24, true, false));
    terms.push_back(Symbol::Symbol("John", 25, true, false));
    terms.push_back(Symbol::Symbol("Mary", 26, true, false));
    terms.push_back(Symbol::Symbol("the", 27, true, false));
    terms.push_back(Symbol::Symbol("man", 28, true, false));
    terms.push_back(Symbol::Symbol("child", 29, true, false));
    terms.push_back(Symbol::Symbol("that", 30, true, false));

    std::vector<std::vector<Symbol::Symbol>> words;
    g.non_terminals.clear();
    g.rules.clear();

    // Não-terminais
    g.non_terminals.push_back(Symbol::Symbol("S", 0, false, false));
    g.non_terminals.push_back(Symbol::Symbol("380", 1, false, false));
    g.non_terminals.push_back(Symbol::Symbol("0", 2, false, false));
    g.non_terminals.push_back(Symbol::Symbol("1", 3, false, false));
    g.non_terminals.push_back(Symbol::Symbol("2", 4, false, false));
    g.non_terminals.push_back(Symbol::Symbol("3", 5, false, false));
    g.non_terminals.push_back(Symbol::Symbol("17", 6, false, false));
    g.non_terminals.push_back(Symbol::Symbol("18", 7, false, false));
    g.non_terminals.push_back(Symbol::Symbol("21", 8, false, false));
    g.non_terminals.push_back(Symbol::Symbol("22", 9, false, false));
    g.non_terminals.push_back(Symbol::Symbol("23", 10, false, false));
    g.non_terminals.push_back(Symbol::Symbol("47", 11, false, false));
    g.non_terminals.push_back(Symbol::Symbol("53", 12, false, false));
    g.non_terminals.push_back(Symbol::Symbol("57", 13, false, false));
    g.non_terminals.push_back(Symbol::Symbol("80", 14, false, false));
    g.non_terminals.push_back(Symbol::Symbol("279", 15, false, false));
    g.non_terminals.push_back(Symbol::Symbol("283", 16, false, false));
    g.non_terminals.push_back(Symbol::Symbol("45", 17, false, false));

    std::vector<Symbol::Symbol> left;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;

    // S -> NP Vi ADV
    left.push_back(g.non_terminals[0]); // S
    rhs.first.push_back(g.non_terminals[1]); // NP
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // NP -> NPa
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[1]); // NP
    rhs.first.push_back(g.non_terminals[5]); // NPa
    rhs.first.push_back(g.non_terminals[13]); // NPa
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // 0 ->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[2]); // Vi
    rhs.first.push_back(g.non_terminals[16]); // 'appears'
    rhs.first.push_back(g.non_terminals[2]); // 'appears'
    rhs.second.first = 0.132338118574;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[5]); // 'appears'
    rhs.first.push_back(g.non_terminals[17]); // 'appears'
    rhs.second.first = 0.0270397939482;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[27]); // 'is'
    rhs.second.first = 0.840622087477;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    //1 ->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[3]); // Vi
    rhs.first.push_back(terms[30]); // 'is'
    rhs.second.first = 1;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    //2->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[4]); // Vi
    rhs.first.push_back(g.non_terminals[15]); // 'appears'
    rhs.first.push_back(g.non_terminals[5]); // 'appears'
    rhs.second.first = 0.000000000000000000000000677392663491;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[2]); // 'seems'
    rhs.second.first = 0.166666666667;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[4]); // 'appears'
    rhs.first.push_back(g.non_terminals[14]); // 'appears'
    rhs.second.first = 0.00000000000000000000000114477054448;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[3]); // 'seems'
    rhs.second.first = 0.166666666667;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[29]); // 'seems'
    rhs.second.first = 0.166666666667;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[28]); // 'seems'
    rhs.second.first = 0.166666666667;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[1]); // 'seems'
    rhs.second.first = 0.166666666667;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[4]); // 'seems'
    rhs.second.first = 0.166666666667;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // 3->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[5]); // ADV
    rhs.first.push_back(g.non_terminals[2]); // 'appears'
    rhs.first.push_back(g.non_terminals[4]); // 'appears'
    rhs.second.first = 0.457527291952;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[5]); // 'appears'
    rhs.first.push_back(g.non_terminals[14]); // 'appears'
    rhs.second.first = 0.122970472532;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[25]); // 'ugly'
    rhs.second.first = 0.209751117758;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[26]); // 'ugly'
    rhs.second.first = 0.209751117758;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // NPa -> 'John' | 'Mary' | 'the' 'man' | 'the' 'child'
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[6]); // NPa
    rhs.first.push_back(terms[7]); // 'John'
    rhs.second.first = 0.25;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[6]); // 'John'
    rhs.second.first = 0.25;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[5]); // 'John'
    rhs.second.first = 0.25;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[8]); // 'John'
    rhs.second.first = 0.25;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // 18 ->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[7]); // VPa
    rhs.first.push_back(g.non_terminals[14]); // Vt
    rhs.first.push_back(g.non_terminals[7]); // NP
    rhs.second.first = 0.0000000110956753618;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[17]); // 'John'
    rhs.second.first = 0.249999997226;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[20]); // 'John'
    rhs.second.first = 0.249999997226;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[18]); // 'John'
    rhs.second.first = 0.249999997226;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[19]); // 'John'
    rhs.second.first = 0.249999997226;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // 21 ->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[8]); // VPa
    rhs.first.push_back(terms[16]); // 'John'
    rhs.second.first = 0.25;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[14]); // 'John'
    rhs.second.first = 0.25;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[13]); // 'John'
    rhs.second.first = 0.25;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[15]); // 'John'
    rhs.second.first = 0.25;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    //22->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[9]); // VPa
    rhs.first.push_back(g.non_terminals[14]); // Vt
    rhs.first.push_back(g.non_terminals[9]); // NP
    rhs.second.first = 0.000000022061171376;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[9]); // Vt
    rhs.first.push_back(g.non_terminals[16]); // NP
    rhs.second.first = 0.00469639274064;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[11]); // 'John'
    rhs.second.first = 0.2488258963;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[10]); // 'John'
    rhs.second.first = 0.2488258963;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[12]); // 'John'
    rhs.second.first = 0.2488258963;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[9]); // 'John'
    rhs.second.first = 0.2488258963;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));


    // 23 ->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[10]); // VPa
    rhs.first.push_back(g.non_terminals[8]); // Vt
    rhs.first.push_back(g.non_terminals[3]); // NP
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));



    // 47 ->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[11]); // VPa
    rhs.first.push_back(terms[21]); // 'John'
    rhs.second.first = 0.333333333333;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[23]); // 'John'
    rhs.second.first = 0.333333333333;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(terms[24]); // 'John'
    rhs.second.first = 0.333333333333;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    //53->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[12]); // VPa
    rhs.first.push_back(terms[22]); // 'John'
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    //57->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[13]); // VPa
    rhs.first.push_back(g.non_terminals[9]); // Vt
    rhs.first.push_back(g.non_terminals[5]); // NP
    rhs.second.first = 0.483;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[7]); // Vt
    rhs.first.push_back(g.non_terminals[6]); // NP
    rhs.second.first = 0.517;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    //80->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[14]); // VPa
    rhs.first.push_back(g.non_terminals[10]); // Vt
    rhs.first.push_back(g.non_terminals[5]); // NP
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    //279->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[15]); // VPa
    rhs.first.push_back(g.non_terminals[4]); // Vt
    rhs.first.push_back(g.non_terminals[12]); // NP
    rhs.second.first = 0.91639350296;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[4]); // Vt
    rhs.first.push_back(g.non_terminals[11]); // NP
    rhs.second.first = 0.0836064970398;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    //279->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[16]); // VPa
    rhs.first.push_back(g.non_terminals[5]); // Vt
    rhs.first.push_back(g.non_terminals[12]); // NP
    rhs.second.first = 0.766393442623;
    right.push_back(rhs);

    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[5]); // Vt
    rhs.first.push_back(g.non_terminals[11]); // NP
    rhs.second.first = 0.233606557377;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // 45 ->
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[17]); // VPa
    rhs.first.push_back(g.non_terminals[10]); // Vt
    rhs.first.push_back(g.non_terminals[2]); // NP
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // Imprime a gramática
    g.terminals = terms;
    g.print_grammar();
    /*std::ofstream outFile("/home/petcomp/CLionProjects/grammatical-inferences/examples/palavrasNL8_" + to_string(max_length) + ".txt");
    words = g.generate_max_size_words_from_rules(max_length);
    cout << words.size() << " " << terms.size() << endl;
    outFile << words.size() << " " << terms.size() << std::endl;
    for (auto w: words) {
        cout << w.size() << " " << g.convert_vector_to_string(w) << endl;
        outFile << w.size() << " " << g.convert_vector_to_string(w) << std::endl;
    }
    outFile.close();
    return words;*/
}
void generate_AC1(Grammar::Grammar &g) {
    std::vector<Symbol::Symbol> terms;
    terms.push_back(Symbol::Symbol("0", 0, true, false));
    terms.push_back(Symbol::Symbol("1", 1, true, false));

    vector<vector<Symbol::Symbol>> words;
    g.non_terminals.clear();
    g.rules.clear();

    g.non_terminals.push_back(Symbol::Symbol("S", 0, false, false));
    g.non_terminals.push_back(Symbol::Symbol("P", 1, false, false));

    vector<Symbol::Symbol> left;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;


    //S -> P
    left.push_back(g.non_terminals[0]);  // S
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // S -> aSbA
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[1]);  // S
    rhs.first.push_back(terms[0]);  // a
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.first.push_back(terms[1]); //b
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);

    //S -> bSaS
    rhs.first.clear();
    rhs.first.push_back(terms[1]);  // a
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.first.push_back(terms[0]); //b
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);

    //S -> lambda
    rhs.first.clear();
    //rhs.first.push_back(terms[0]);
    rhs.first.push_back(Symbol::Symbol("", -1, true, false));  // lambda
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    /*rhs.first.clear();
    rhs.first.push_back(terms[1]);
    //rhs.first.push_back(Symbol::Symbol("", -1, true, false));  // lambda
    rhs.second.first = 1.0/4.0;
    right.push_back(rhs);*/


    g.terminals = terms;

}


void generate_AC2(Grammar::Grammar &g) {
    std::vector<Symbol::Symbol> terms;
    terms.push_back(Symbol::Symbol("0", 0, true, false));
    terms.push_back(Symbol::Symbol("1", 1, true, false));

    vector<vector<Symbol::Symbol>> words;
    g.non_terminals.clear();
    g.rules.clear();

    g.non_terminals.push_back(Symbol::Symbol("P", 0, false, false));
    g.non_terminals.push_back(Symbol::Symbol("S", 1, false, false));
    g.non_terminals.push_back(Symbol::Symbol("A", 2, false, false));
    g.non_terminals.push_back(Symbol::Symbol("B", 3, false, false));
    g.non_terminals.push_back(Symbol::Symbol("C", 4, false, false));

    vector<Symbol::Symbol> left;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;


    //P -> S
    left.push_back(g.non_terminals[0]);  // S
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.second.first = 1.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // S -> BC
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[1]);  // S
    rhs.first.push_back(g.non_terminals[3]);  // S
    rhs.first.push_back(g.non_terminals[4]);  // S
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);

    //S -> CB
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[4]);  // S
    rhs.first.push_back(g.non_terminals[3]);  // S
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);

    //S -> lambda
    /*rhs.first.clear();
    //rhs.first.push_back(terms[0]);
    rhs.first.push_back(Symbol::Symbol("", -1, true, false));  // lambda
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);*/
    g.rules.push_back(Rule::Rule(left, right));




    // A -> AS
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[2]);  // S
    rhs.first.push_back(g.non_terminals[2]);  // S
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);

    //A -> SA
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.first.push_back(g.non_terminals[2]);  // S
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);

    //A -> 'a'
    rhs.first.clear();
    rhs.first.push_back(terms[0]);  // a
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // B -> BS
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[3]);  // S
    rhs.first.push_back(g.non_terminals[3]);  // S
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);

    //B -> SB
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[1]);  // S
    rhs.first.push_back(g.non_terminals[3]);  // S
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);

    //B -> 'b'
    rhs.first.clear();
    rhs.first.push_back(terms[1]);  // b
    rhs.second.first = 1.0/3.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));

    // C -> BA
    left.clear();
    right.clear();
    rhs.first.clear();
    left.push_back(g.non_terminals[4]);  // S
    rhs.first.push_back(g.non_terminals[3]);  // S
    rhs.first.push_back(g.non_terminals[2]);  // S
    rhs.second.first = 1.0/2.0;
    right.push_back(rhs);

    //C -> AB
    rhs.first.clear();
    rhs.first.push_back(g.non_terminals[2]);  // S
    rhs.first.push_back(g.non_terminals[3]);  // S
    rhs.second.first = 1.0/2.0;
    right.push_back(rhs);
    g.rules.push_back(Rule::Rule(left, right));


    g.n_non_terminals = g.non_terminals.size();
    g.terminals = terms;

}
