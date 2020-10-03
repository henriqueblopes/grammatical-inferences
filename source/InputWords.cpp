//
// Created by henrique on 04/05/2020.
//

#include "InputWords.h"
#include <string>
#include <iostream>
#include <sstream>
#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
namespace fs = std::experimental::filesystem;

InputWords::InputWords(bool timed, int nTerminals) : timed(timed), nTerminals(nTerminals) {}

void InputWords::readWords() {
    fs::path p = fs::current_path().parent_path();
    fs::path p1 = (p /= "McGill-Billboard");
    for(auto& p: fs::recursive_directory_iterator(p1)) {
        if (!fs::is_directory(p)) {
            //std::cout << p.path() << '\n';
            convertFiletoWord(p.path(), 7);
        }
    }
    //cout << "nSequences: " << inputWords.size() << endl;
}

//TO DO Ver como Tratar o acorde N
void InputWords::convertFiletoWord(fs::path path, int minSize) {
    string targeTone = "C";
    std::ifstream ifs (path, std::ifstream::in);
    string s;
    getline(ifs,s);
    vector<Symbol> timedWord;
    vector<Symbol> noTimedWord;
    vector<Symbol> timedChordLine;
    vector<Symbol> noTimedChordLine;
    int beat = 0 ;
    string tone = "";
    while (ifs.good()) {
        //std::cout << s << endl;
        size_t tokenPos = s.find("metre: ");
        if (tokenPos != string::npos) {
            beat = stoi(s.substr(tokenPos+7, s.find("/") - (tokenPos+7)));
            //cout << "nBeats:" << beat << endl;
        }
        tokenPos = s.find("tonic: ");
        if (tokenPos != string::npos) {
            if (!noTimedWord.empty()) {
                if (noTimedWord.size() >= minSize) {
                    transposeTo(tone, targeTone, noTimedWord);
                    transposeTo(tone, targeTone, timedWord);
                    if(timed)
                        inputWords.push_back(timedWord);
                    else inputWords.push_back(noTimedWord);
                }
                noTimedWord.clear();
                timedWord.clear();
            }
            tone = s.substr(tokenPos+7, s.size() - (tokenPos+7));
//            cout << "Tone:" << tone << endl;
        }
        if (!tone.empty()) {
            size_t initPipe = s.find("|");
            if (initPipe != string::npos) {
                size_t comma = s.find(",");
                comma = s.find(",",comma+1);
                if (comma != string::npos) {
                    if (comma < initPipe) {
                        if (!noTimedWord.empty()) {
                            if (noTimedWord.size() >= minSize) {
                                transposeTo(tone, targeTone, noTimedWord);
                                transposeTo(tone, targeTone, timedWord);
                                if(timed)
                                    inputWords.push_back(timedWord);
                                else inputWords.push_back(noTimedWord);
                            }
                            noTimedWord.clear();
                            timedWord.clear();
                        }

                    }
                }
                while (initPipe != string::npos) {
                    size_t nextPipe = s.find("|", initPipe+1);
                    if (nextPipe == string::npos && initPipe != string::npos && initPipe+1 < s.size()) {
                        if (s.substr(initPipe+2,initPipe+3 - (initPipe+2)).compare("x") == 0) {
                            int times = stoi(s.substr(initPipe+3, 1));
                            vector<Symbol> aux = noTimedChordLine;
                            vector<Symbol> aux2 = timedChordLine;
                            for (int i = 0; i < times; i++) {
                                noTimedChordLine.insert(noTimedChordLine.end(), aux.begin(), aux.end());
                                timedChordLine.insert(timedChordLine.end(), aux2.begin(), aux2.end());
                            }
                        }
                    }

                    size_t initBlank = s.find(" ", initPipe+1);
                    int chordAmanout = 0;
                    int tempBeat = beat;
                    vector<Symbol> chords;
                    vector<Symbol> timedChords;
                    string currentSymbol = "";
                    string currentChord = "";
                    if(initBlank != string::npos && nextPipe != string::npos) {
                        while (initBlank+1 < nextPipe) {
                            size_t nextBlank = s.find(" ", initBlank+1);
                            currentSymbol = s.substr(initBlank+1, nextBlank - (initBlank+1));
                            if (currentSymbol.compare(".") == 0) {
                                chordAmanout++;
                                timedChords.push_back(Symbol(currentChord, 0, true, false));
                            }
                            else if(currentSymbol.substr(0,1).compare("(") == 0)
                                tempBeat = stoi(s.substr(initPipe+3, s.find("/") - (initPipe+3)));
                            else if (currentSymbol.compare("|") != 0){
                                chordAmanout++;
                                currentChord = currentSymbol;
                                if (currentSymbol.compare("") == 0)
                                    cout << "some error" << endl;
                                chords.push_back(Symbol(currentSymbol, 0, true, false));
                                timedChords.push_back(Symbol(currentSymbol, 0, true, false));
                            }
                            initBlank = nextBlank;
                        }
                        vector<Symbol>::iterator itTimedChords;
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
                    transposeTo(tone, targeTone, noTimedWord);
                    transposeTo(tone, targeTone, timedWord);
                    if(timed)
                        inputWords.push_back(timedWord);
                    else inputWords.push_back(noTimedWord);
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

void InputWords::transposeTo(string actualTone, string targetTone, vector<Symbol> & word) {
    std::vector<string> tones = {"C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"};
    std::vector<string> tonesS = {"C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"};
    while (tones[0].compare(targetTone) != 0)
        rotate(tones.begin(), tones.begin()+1, tones.end());
    int posTone;
    for (int i = 0; i< tones.size(); i++) {
        if (tones[i].compare(actualTone) == 0) {
            posTone = tones.size() - i;
            break;
        } else if (tonesS[i].compare(actualTone) == 0) {
            posTone = tones.size() -i;
            break;
        }
    }
    for ( int j = 0; j < word.size(); j++) {
        for (int i = 0; i < tones.size(); i ++) {
            if (word[j].name.substr(1,1).compare("b") != 0  && word[j].name.substr(1,1).compare("#") != 0) {
                if (word[j].name.substr(0,1).compare(tones[i]) == 0)     {
                    word[j].name = tones[(i+posTone)%tones.size()] + word[j].name.substr(1,word[j].name.size()-1);
                    break;
                }
            }
            else if (word[j].name.substr(0,2).compare(tones[i]) == 0 || word[j].name.substr(0,2).compare(tonesS[i]) == 0) {
                string aux = tones[(i+posTone)%tones.size()] + word[j].name.substr(2,word[j].name.size()-2);
                word[j].name = tones[(i+posTone)%tones.size()] + word[j].name.substr(2,word[j].name.size()-2);
                break;
            }
        }

    }
}

void InputWords::countChords() {
    for (vector<Symbol> w: inputWords) {
        for (Symbol s: w) {
            if (chordCounts.find(s.name) == chordCounts.end())
                chordCounts[s.name] = 1;
            else
                chordCounts[s.name]++;
        }
    }
}

vector<Symbol> InputWords::generateTerminals() {
    countChords();
    vector<Symbol> terminals;
    pair<string,int> maxChord;
    int i =0;
    while (i < nTerminals-1) {
        maxChord = make_pair("", 0);
        for (auto a: chordCounts) {
            if(a.second > maxChord.second)
                maxChord = a;
        }
        Symbol aux = Symbol(maxChord.first, i, true, false);
        terminals.push_back(aux);
        chordMap[aux.name] = aux;
        chordCounts.erase(maxChord.first);
        i++;
    }
    Symbol aux = Symbol("Other", nTerminals-1, true, false);
    terminals.push_back(aux);
    chordMap["Other"] = aux;
    truncChordWords();
    return terminals;
}

void InputWords::truncChordWords() {
    for (vector<vector<Symbol>>::iterator  itIW = inputWords.begin(); itIW != inputWords.end(); itIW++) {
        for (vector<Symbol>::iterator itS = (*itIW).begin(); itS != (*itIW).end(); itS++) {
            if (chordMap.find((*itS).name) == chordMap.end()) {
                (*itS) = chordMap["Other"];
            } else {
                (*itS) = chordMap[(*itS)
                                  .name];
            }
        }
    }
}



void InputWords::checkEmptyString(string path) {
    for (auto w: inputWords) {
        for (auto s: w) {
            if (s.name.compare("") ==0) {
                cout <<" some error in " << path << endl;
                break;
            }
        }
    }
}

void InputWords::printWordSizes() {
    double sum = 0.0;
    int max = 0;
    for (auto w: inputWords) {
        //cout << "Size:  " << w.size() << endl;
        sum += w.size();
        if (w.size() > max)
            max = w.size();
    }
    cout << "MaxSize: " << max << " MeanSize: " << sum/inputWords.size() << endl;
}

void InputWords::selectTrainingWords(int nSharesOrAmount, bool byShare) {
    if (byShare)
        nTestShares = nSharesOrAmount;
    else
        nTestShares = inputWords.size()/nSharesOrAmount;
    actualShare = 0;
    //std::random_shuffle ( inputWords.begin(), inputWords.end());
    int tAmount = (int) inputWords.size()/nTestShares;
    vector<vector<Symbol>> tWords;
    for (int i = 0; i< tAmount; i++) {
        testWords.push_back(inputWords[0]);
        inputWords.erase(inputWords.begin());
    }
}

bool InputWords::nextShareTrainingWords() {
    actualShare++;
    if (actualShare >= nTestShares)
        return false;
    inputWords.insert(inputWords.end(), testWords.begin(), testWords.end());
    testWords.clear();
    int tAmount = (int) inputWords.size()/nTestShares;
    for (int i = 0; i< tAmount; i++) {
        testWords.push_back(inputWords[0]);
        inputWords.erase(inputWords.begin());
    }
    return true;
}

void InputWords::iterateChords() {
    string minor = "10110101101";
    string major = "10101101010";
    string chordNotes = "00000000000";
    countChords();
    std::unordered_map<string, int>::iterator chordCountsIt;
    size_t tokenPos;
    string tone = "";
    string mode = "";
    string bass = "";
    string addition = "";
    for (chordCountsIt = chordCounts.begin(); chordCountsIt != chordCounts.end(); chordCountsIt++) {
        //cout << (*chordCountsIt).first << endl;
        if ((*chordCountsIt).first.compare("N") !=0 && (*chordCountsIt).first.compare("&pause") !=0 && (*chordCountsIt).first.compare("*") !=0) {
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


            chordNotes = buildChordVector(tone, mode, addition, bass);

            std::vector<string> tones = {"C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"};
            std::vector<string> tonesS = {"B#", "C#", "D", "D#", "Fb", "E#", "F#", "G", "G#", "A", "A#", "Cb"};
            while (tones[0].compare(tone) != 0 && tonesS[0].compare(tone) != 0)  {
                rotate(chordNotes.begin(), chordNotes.begin()+1, chordNotes.end());
                rotate(tones.begin(), tones.begin()+1, tones.end());
                rotate(tonesS.begin(), tonesS.begin()+1, tonesS.end());
            }

            chordMapReduction.insert(make_pair((*chordCountsIt).first, chordNotes));
        } else {
            chordMapReduction.insert(make_pair((*chordCountsIt).first, "00000000000"));
        }
    }
    countAndShowReductedChords();

}
string InputWords::buildChordVector(string tone, string mode, string addition, string bass) {
    string chord = "00000000000";
    chord[0] = '1';

    if (mode.size() > 2) {    //ModeSize > 2
        if (mode.substr(0, 3).compare("maj") == 0) {
            chord[4] = chord[7] = '1';
        } else if (mode.substr(0, 3).compare("min") == 0) {
            chord[3] = chord[7] = '1';
        } else if (mode.substr(0, 3).compare("dim") == 0) {
            chord[3] = chord[6] = '1';
        } else if (mode.substr(0, 3).compare("aug") == 0) {
            chord[4] = chord[8] = '1';
        }
        if (mode.size() >= 4) {
            if (mode.substr(mode.size() - 4, 4).compare("maj7") == 0) {
                chord[11] = '1';
            } else if (mode.substr(mode.size() - 4, mode.size()).compare("min7") == 0) {
                chord[10] = '1';
            } else if (mode.compare("dim7") == 0) {
                chord[9] = '1';
            } else if (mode.compare("hdim7") == 0) {
                chord[3] = chord[6] = chord[10] = '1';
            } else if (mode.substr(mode.size() - 1, mode.size()).compare("6") == 0) {
                chord[9] = '1';
            } else if (mode.compare("maj9") == 0) {
                chord[11] = chord[2] = '1';
            } else if (mode.compare("min9") == 0) {
                chord[10] = chord[2] = '1';
            }
                //Talvez omitir notas no 11 e 13
            else if (mode.compare("maj11") == 0) {
                chord[11] = chord[2] = chord[5] = '1';
            } else if (mode.compare("min11") == 0) {
                chord[10] = chord[2] = chord[5] = '1';
            } else if (mode.compare("maj13") == 0) {
                chord[11] = chord[2] = chord[5] = chord[9] = '1';
            } else if (mode.compare("min13") == 0) {
                chord[10] = chord[2] = chord[5] = chord[9] = '1';
            } else if (mode.compare("sus4") == 0) {
                chord[6] = chord[7] = '1';
            } else if (mode.compare("sus2") == 0) {
                chord[2] = chord[7] = '1';
            }
        }

    } else {
        if (mode.compare("5") == 0) {
            chord[7] = '1';
        } else {
            chord[4] = chord[7] = '1';
            if (mode.compare ("7") == 0) {
                chord[10] = '1';
            } else if (mode.compare ("9") == 0) {
                chord[10] = chord[2] = '1';
            } else if (mode.compare ("11") == 0) {
                chord[10] = chord[2] = chord[5] = '1';
            } else if (mode.compare ("13") == 0) {
                chord[10] = chord[2] = chord[5] = chord[9] = '1';
            }
        }
    }

    if (addition.size()>0) {
        size_t commaPos = addition.find(",");
        string note;
        size_t commaBefore = 0;
        do {
            if (commaPos == string::npos)
                note = addition.substr(commaBefore, addition.size()-commaBefore);
            else
                note = addition.substr(commaBefore, commaPos-commaBefore);
            chord = addChordNote(chord, note);
            commaBefore = commaPos+1;
            commaPos = addition.find(",", commaPos+1);
        } while (commaPos != string::npos);
    }

    if(bass.size() >0) {
        chord = addChordNote(chord, bass);
    }
    return chord;
}

string InputWords::addChordNote(string chord, string note) {
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

void InputWords::countAndShowReductedChords() {
    std::unordered_map<string, string>::iterator chordMRIt;
    for (chordMRIt = chordMapReduction.begin(); chordMRIt != chordMapReduction.end(); chordMRIt++) {
        if (reductedChordCounts.find((*chordMRIt).second) == chordCounts.end())
            reductedChordCounts[(*chordMRIt).second] = 1;
        else
            reductedChordCounts[(*chordMRIt).second]++;
    }
    for (auto a: reductedChordCounts) {
        cout <<a.first << " " << a.second << endl;
    }
}
