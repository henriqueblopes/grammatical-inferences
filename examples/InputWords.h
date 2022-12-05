//
// Created by henrique on 04/05/2020.
//

#ifndef GRAMMARINDUCTION_INPUTWORDS_H
#define GRAMMARINDUCTION_INPUTWORDS_H


#include <vector>
#include <gInfer/Symbol.h>
#ifdef CXX_FILESYSTEM_IS_EXPERIMENTAL
    #include <experimental/filesystem>
#else
    #include <filesystem>
#endif
#include <iostream>
#include <fstream>
#include <unordered_map>

#ifdef CXX_FILESYSTEM_IS_EXPERIMENTAL
    namespace fs = std::experimental::filesystem;
#else
    namespace fs = std::filesystem;
#endif


using namespace std;

class InputWords {
public:
    vector<vector<Symbol::Symbol>> input_words;
    vector<vector<Symbol::Symbol>> test_words;
    int actual_share;
    std::unordered_map<string, int> reducted_chord_counts;
    std::unordered_map<string, int> chord_counts;

private:
    bool timed;

    std::unordered_map<string, Symbol::Symbol> chord_map;
    std::unordered_map<string, string> chord_map_reduction;

    int n_terminals;
    int n_test_shares;

public:
    InputWords(bool timed, int n);
    void read_words(bool is_entire_music);
    void convert_file_to_word (const fs::path& path, unsigned long minSize);
    static void transpose_to(const string& actual_tone, const string& target_tone, vector<Symbol::Symbol> & word);
    vector<Symbol::Symbol> generate_terminals(std::unordered_map<string, int> counted_chords);
    void print_word_sizes();
    void select_training_words(int n_shares_or_amount, bool by_share);
    void count_chords();
    void trunc_chord_words();
    void check_empty_string(const string& path);
    bool next_share_training_words();
    void iterate_chords();
    static string build_chord_vector(const string& mode, const string& addition, const string& bass);
    void change_words_to_reducted_chords();
private:
    static string add_chord_note(string chord, string note);

    void count_and_show_reducted_chords();
};



#endif //GRAMMARINDUCTION_INPUTWORDS_H
