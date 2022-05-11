//
// Created by henrique on 04/05/2020.
//

#include "InputWords.h"
#include <string>
#include <algorithm>


InputWords::InputWords(bool timed, int n) : timed(timed), n_terminals(n) {}

void InputWords::read_words() {
    fs::path p = fs::current_path().parent_path();
    fs::path p1 = (p /= "resources/McGill-Billboard");
    for(auto& pl: fs::recursive_directory_iterator(p1)) {
        if (!fs::is_directory(pl)) {
            //std::cout << p.path() << '\n';
            convert_file_to_word(pl.path(), 7);
        }
    }
    //cout << "nSequences: " << input_words.size() << endl;
}

//TODO Ver como Tratar o acorde N
void InputWords::convert_file_to_word(const fs::path& path, unsigned long minSize) {
    string targeTone = "C";
    std::ifstream ifs (path, std::ifstream::in);
    string s;
    getline(ifs,s);
    vector<Symbol::Symbol> timedWord;
    vector<Symbol::Symbol> noTimedWord;
    vector<Symbol::Symbol> timedChordLine;
    vector<Symbol::Symbol> noTimedChordLine;
    int beat = 0 ;
    string tone;
    tone = "";
    while (ifs.good()) {
        //std::cout << s << endl;
        size_t tokenPos = s.find("metre: ");
        if (tokenPos != string::npos) {
            beat = stoi(s.substr(tokenPos+7, s.find('/') - (tokenPos+7)));
            //cout << "nBeats:" << beat << endl;
        }
        tokenPos = s.find("tonic: ");
        if (tokenPos != string::npos) {
            if (!noTimedWord.empty()) {
                if (noTimedWord.size() >= minSize) {
                    transpose_to(tone, targeTone, noTimedWord);
                    transpose_to(tone, targeTone, timedWord);
                    if(timed)
                        input_words.push_back(timedWord);
                    else input_words.push_back(noTimedWord);
                }
                noTimedWord.clear();
                timedWord.clear();
            }
            tone = s.substr(tokenPos+7, s.size() - (tokenPos+7));
//            cout << "Tone:" << tone << endl;
        }
        if (!tone.empty()) {
            size_t initPipe = s.find('|');
            if (initPipe != string::npos) {
                size_t comma = s.find(',');
                comma = s.find(',',comma+1);
                if (comma != string::npos) {
                    if (comma < initPipe) {
                        if (!noTimedWord.empty()) {
                            if (noTimedWord.size() >= minSize) {
                                transpose_to(tone, targeTone, noTimedWord);
                                transpose_to(tone, targeTone, timedWord);
                                if(timed)
                                    input_words.push_back(timedWord);
                                else input_words.push_back(noTimedWord);
                            }
                            noTimedWord.clear();
                            timedWord.clear();
                        }

                    }
                }
                while (initPipe != string::npos) {
                    size_t nextPipe = s.find('|', initPipe+1);
                    if (nextPipe == string::npos && initPipe != string::npos && initPipe+1 < s.size()) {
                        if (s.substr(initPipe+2,initPipe+3 - (initPipe+2)) == "x") {
                            int times = stoi(s.substr(initPipe+3, 1));
                            vector<Symbol::Symbol> aux = noTimedChordLine;
                            vector<Symbol::Symbol> aux2 = timedChordLine;
                            for (int i = 0; i < times; i++) {
                                noTimedChordLine.insert(noTimedChordLine.end(), aux.begin(), aux.end());
                                timedChordLine.insert(timedChordLine.end(), aux2.begin(), aux2.end());
                            }
                        }
                    }

                    size_t initBlank = s.find(' ', initPipe+1);
                    int chordAmanout = 0;
                    int tempBeat = beat;
                    vector<Symbol::Symbol> chords;
                    vector<Symbol::Symbol> timedChords;
                    string currentSymbol;
                    string currentChord;
                    if(initBlank != string::npos && nextPipe != string::npos) {
                        while (initBlank+1 < nextPipe) {
                            size_t nextBlank = s.find(' ', initBlank+1);
                            currentSymbol = s.substr(initBlank+1, nextBlank - (initBlank+1));
                            if (currentSymbol == ".") {
                                chordAmanout++;
                                timedChords.emplace_back(currentChord, 0, true, false);
                            }
                            else if(currentSymbol.substr(0,1) == "(")
                                tempBeat = stoi(s.substr(initPipe+3, s.find('/') - (initPipe+3)));
                            else if (currentSymbol != "|"){
                                chordAmanout++;
                                currentChord = currentSymbol;
                                if (currentSymbol.empty())
                                    cout << "some error" << endl;
                                chords.emplace_back(currentSymbol, 0, true, false);
                                timedChords.emplace_back(currentSymbol, 0, true, false);
                            }
                            initBlank = nextBlank;
                        }
                        vector<Symbol::Symbol>::iterator itTimedChords;
                        for (itTimedChords = timedChords.begin(); itTimedChords != timedChords.end(); itTimedChords++)
                            for (int i = 0; i < tempBeat/chordAmanout; i++)
                                timedChordLine.push_back((*itTimedChords));
                        noTimedChordLine.insert(noTimedChordLine.end(), chords.begin(), chords.end());
                    }
                    initPipe = nextPipe;
                }
            }
        }
        if (s.find("end") != string::npos) {
            if (!noTimedWord.empty()) {
                if (noTimedWord.size() >= minSize) {
                    transpose_to(tone, targeTone, noTimedWord);
                    transpose_to(tone, targeTone, timedWord);
                    if(timed)
                        input_words.push_back(timedWord);
                    else input_words.push_back(noTimedWord);
                }
                noTimedWord.clear();
                timedWord.clear();
            }
        }
        noTimedWord.insert(noTimedWord.end(), noTimedChordLine.begin(), noTimedChordLine.end());
        timedWord.insert(timedWord.end(), timedChordLine.begin(), timedChordLine.end());
        noTimedChordLine.clear();
        timedChordLine.clear();
        getline(ifs,s);
    }
    ifs.close();
}

void InputWords::transpose_to(const string& actual_tone, const string& target_tone, vector<Symbol::Symbol> & word) {
    std::vector<string> tones = {"C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"};
    std::vector<string> tonesS = {"C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"};
    while (tones[0] != target_tone)
        rotate(tones.begin(), tones.begin()+1, tones.end());
    int posTone = 0;
    for (unsigned long i = 0; i< tones.size(); i++) {
        if (tones[i] == actual_tone || tonesS[i] == actual_tone) {
            posTone = (int) (tones.size() - i);
            break;
        }
    }
    for (auto & j : word) {
        for (unsigned long i = 0; i < tones.size(); i ++) {
            if (j.name.substr(1,1) != "b"  && j.name.substr(1,1) != "#") {
                if (j.name.substr(0,1) == tones[i])     {
                    j.name = tones[(i+posTone)%tones.size()] + j.name.substr(1,j.name.size()-1);
                    break;
                }
            }
            else if (j.name.substr(0,2) == tones[i] || j.name.substr(0,2) == tonesS[i]) {
                string aux = tones[(i+posTone)%tones.size()] + j.name.substr(2,j.name.size()-2);
                j.name = tones[(i+posTone)%tones.size()] + j.name.substr(2,j.name.size()-2);
                break;
            }
        }

    }
}

void InputWords::count_chords() {
    for (const vector<Symbol::Symbol>& w: input_words) {
        for (const Symbol::Symbol& s: w) {
            if (chord_counts.find(s.name) == chord_counts.end())
                chord_counts[s.name] = 1;
            else
                chord_counts[s.name]++;
        }
    }
}

vector<Symbol::Symbol> InputWords::generate_terminals(std::unordered_map<string, int> counted_chords) {
    //count_chords();
    vector<Symbol::Symbol> terminals;
    pair<string,int> maxChord;
    int i =0;
    while (i < n_terminals - 1) {
        maxChord = make_pair("", 0);
        for (const auto& a: counted_chords) {
            if(a.second > maxChord.second)
                maxChord = a;
        }
        Symbol::Symbol aux = Symbol::Symbol(maxChord.first, i, true, false);
        terminals.push_back(aux);
        chord_map[aux.name] = aux;
        counted_chords.erase(maxChord.first);
        i++;
    }
    Symbol::Symbol aux = Symbol::Symbol("Other", n_terminals - 1, true, false);
    terminals.push_back(aux);
    chord_map["Other"] = aux;
    trunc_chord_words();
    return terminals;
}

void InputWords::trunc_chord_words() {
    for (auto & inputWord : input_words) {
        for (auto & itS : inputWord) {
            if (chord_map.find(itS.name) == chord_map.end()) {
                itS = chord_map["Other"];
            } else {
                itS = chord_map[itS.name];
            }
        }
    }
}



void InputWords::check_empty_string(const string& path) {
    for (const auto& w: input_words) {
        for (const auto& s: w) {
            if (s.name.empty()) {
                cout <<" some error in " << path << endl;
                break;
            }
        }
    }
}

void InputWords::print_word_sizes() {
    double sum = 0.0;
    size_t max = 0;
    for (const auto& w: input_words) {
        //cout << "Size:  " << w.size() << endl;
        sum += w.size();
        if (w.size() > max)
            max = w.size();
    }
    cout << "MaxSize: " << max << " MeanSize: " << sum / input_words.size() << endl;
}

void InputWords::select_training_words(int n_shares_or_amount, bool by_share) {
    if (by_share)
        n_test_shares = n_shares_or_amount;
    else
        n_test_shares = ((int) input_words.size()) / n_shares_or_amount;
    actual_share = 0;
    //std::random_shuffle ( input_words.begin(), input_words.end());
    int tAmount = (int) input_words.size() / n_test_shares;
    vector<vector<Symbol::Symbol>> tWords;
    for (int i = 0; i< tAmount; i++) {
        test_words.push_back(input_words[0]);
        input_words.erase(input_words.begin());
    }
}

bool InputWords::next_share_training_words() {
    actual_share++;
    if (actual_share >= n_test_shares)
        return false;
    input_words.insert(input_words.end(), test_words.begin(), test_words.end());
    test_words.clear();
    int tAmount = (int) input_words.size() / n_test_shares;
    for (int i = 0; i< tAmount; i++) {
        test_words.push_back(input_words[0]);
        input_words.erase(input_words.begin());
    }
    return true;
}

void InputWords::iterate_chords() {
    /*string minor = "10110101101";
    string major = "10101101010";*/
    string chordNotes = "00000000000";
    count_chords();
    std::unordered_map<string, int>::iterator chordCountsIt;
    size_t tokenPos;
    string tone;
    string mode;
    string bass;
    string addition;
    for (chordCountsIt = chord_counts.begin(); chordCountsIt != chord_counts.end(); chordCountsIt++) {
        //cout << (*chordCountsIt).first << endl;
        if ((*chordCountsIt).first !="N" && (*chordCountsIt).first !="&pause" && (*chordCountsIt).first !="*") {
            tokenPos = (*chordCountsIt).first.find(":");
            tone = (*chordCountsIt).first.substr(0, tokenPos);
            tokenPos = (*chordCountsIt).first.find("(");
            if (tokenPos == string::npos) {
                tokenPos = (*chordCountsIt).first.find("/");
                if (tokenPos == string::npos) {
                    mode = (*chordCountsIt).first.substr(tone.size()+1, (*chordCountsIt).first.size()-tone.size()-1);
                }
                else {
                    mode = (*chordCountsIt).first.substr(tone.size()+1, tokenPos-tone.size()-1);
                    bass = (*chordCountsIt).first.substr(tokenPos+1, (*chordCountsIt).first.size()-tokenPos-1);
                }
            } else {
                size_t tokenPosEnd = (*chordCountsIt).first.find(")");
                mode = (*chordCountsIt).first.substr(tone.size()+1, tokenPos-tone.size()-1);
                addition = (*chordCountsIt).first.substr(tokenPos+1, tokenPosEnd-tokenPos-1);
                tokenPos = (*chordCountsIt).first.find("/");
                if (tokenPos != string::npos) {
                    bass = (*chordCountsIt).first.substr(tokenPos+1, (*chordCountsIt).first.size()-tokenPos-1);
                }
            }


            chordNotes = build_chord_vector(mode, addition, bass);

            std::vector<string> tones = {"C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"};
            std::vector<string> tonesS = {"B#", "C#", "D", "D#", "Fb", "E#", "F#", "G", "G#", "A", "A#", "Cb"};
            while (tones[0] != tone && tonesS[0] != tone)  {
                rotate(chordNotes.begin(), chordNotes.begin()+1, chordNotes.end());
                rotate(tones.begin(), tones.begin()+1, tones.end());
                rotate(tonesS.begin(), tonesS.begin()+1, tonesS.end());
            }

            chord_map_reduction.insert(make_pair((*chordCountsIt).first, chordNotes));
        } else {
            chord_map_reduction.insert(make_pair((*chordCountsIt).first, "00000000000"));
        }
    }
    count_and_show_reducted_chords();

}
string InputWords::build_chord_vector(const string& mode, const string& addition, const string& bass) {
    string chord = "00000000000";
    chord[0] = '1';

    if (mode.size() > 2) {    //ModeSize > 2
        if (mode.substr(0, 3) == "maj") {
            chord[4] = chord[7] = '1';
        } else if (mode.substr(0, 3) == "min") {
            chord[3] = chord[7] = '1';
        } else if (mode.substr(0, 3) == "dim") {
            chord[3] = chord[6] = '1';
        } else if (mode.substr(0, 3) == "aug") {
            chord[4] = chord[8] = '1';
        }
        if (mode.size() >= 4) {
            if (mode.substr(mode.size() - 4, 4) == "maj7") {
                chord[11] = '1';
            } else if (mode.substr(mode.size() - 4, mode.size()) == "min7") {
                chord[10] = '1';
            } else if (mode == "dim7" || mode.substr(mode.size() - 1, mode.size()) == "6") {
                chord[9] = '1';
            } else if (mode == "hdim7") {
                chord[3] = chord[6] = chord[10] = '1';
            } else if (mode == "maj9") {
                chord[11] = chord[2] = '1';
            } else if (mode == "min9") {
                chord[10] = chord[2] = '1';
            }
                //Talvez omitir notas no 11 e 13
            else if (mode == "maj11") {
                chord[11] = chord[2] = chord[5] = '1';
            } else if (mode == "min11") {
                chord[10] = chord[2] = chord[5] = '1';
            } else if (mode == "maj13") {
                chord[11] = chord[2] = chord[5] = chord[9] = '1';
            } else if (mode == "min13") {
                chord[10] = chord[2] = chord[5] = chord[9] = '1';
            } else if (mode == "sus4") {
                chord[6] = chord[7] = '1';
            } else if (mode == "sus2") {
                chord[2] = chord[7] = '1';
            }
        }

    } else {
        if (mode == "5") {
            chord[7] = '1';
        } else {
            chord[4] = chord[7] = '1';
            if (mode == "7") {
                chord[10] = '1';
            } else if (mode == "9") {
                chord[10] = chord[2] = '1';
            } else if (mode == "11") {
                chord[10] = chord[2] = chord[5] = '1';
            } else if (mode == "13") {
                chord[10] = chord[2] = chord[5] = chord[9] = '1';
            }
        }
    }

    if (!addition.empty()) {
        size_t commaPos = addition.find(',');
        string note;
        size_t commaBefore = 0;
        do {
            if (commaPos == string::npos)
                note = addition.substr(commaBefore, addition.size()-commaBefore);
            else
                note = addition.substr(commaBefore, commaPos-commaBefore);
            chord = add_chord_note(chord, note);
            commaBefore = commaPos+1;
            commaPos = addition.find(',', commaPos+1);
        } while (commaPos != string::npos);
    }

    if(!bass.empty()) {
        chord = add_chord_note(chord, bass);
    }
    return chord;
}

string InputWords::add_chord_note(string chord, string note) {
    vector<int> intervals = {0,2,4,5,7,9,11};
    int accident = 0;
    while (note[0] =='b') {
        accident -=1;
        note = note.substr(1,note.size()-1);
    }
    while (note[0] =='#') {
        accident +=1;
        note = note.substr(1,note.size()-1);
    }

    int degree = stoi(note);
    chord[intervals[((degree+accident)%7)-1]] = '1';
    return chord;
}

void InputWords::count_and_show_reducted_chords() {
    std::unordered_map<string, string>::iterator chordMRIt;
    for (chordMRIt = chord_map_reduction.begin(); chordMRIt != chord_map_reduction.end(); chordMRIt++) {
        if (reducted_chord_counts.find((*chordMRIt).second) == chord_counts.end())
            reducted_chord_counts[(*chordMRIt).second] = 1;
        else
            reducted_chord_counts[(*chordMRIt).second]++;
    }
    for (const auto& a: reducted_chord_counts) {
        cout <<a.first << " " << a.second << endl;
    }
}
void InputWords::change_words_to_reducted_chords() {
    for (auto  & w: input_words)
        for (auto & s: w)
            s.name = chord_map_reduction[s.name];
}
