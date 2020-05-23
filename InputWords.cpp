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
    fs::path p1 = "/home/henrique/CLionProjects/GrammarInduction/McGill-Billboard"; // portable format
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
        terminals.push_back(Symbol(maxChord.first, i, true, false));
        chordCounts.erase(maxChord.first);
        i++;
    }
    terminals.push_back(Symbol("Other", nTerminals-1, true, false));
    truncChordWords();
    return terminals;
}

void InputWords::truncChordWords() {
    for (vector<vector<Symbol>>::iterator  itIW = inputWords.begin(); itIW != inputWords.end(); itIW++) {
        for (vector<Symbol>::iterator itS = (*itIW).begin(); itS != (*itIW).end(); itS++) {
            if (chordCounts.find((*itS).name) != chordCounts.end()) {
                (*itS).name = "Other";
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
    std::random_shuffle ( inputWords.begin(), inputWords.end());
    int tAmount = (int) inputWords.size()/nTestShares;
    vector<vector<Symbol>> tWords;
    for (int i = 0; i< tAmount; i++) {
        testWords.push_back(inputWords[0]);
        inputWords.erase(inputWords.begin());
    }
}

bool InputWords::nextShareTrainingWords() {
    if (actualShare >= nTestShares)
        return false;
    actualShare++;
    inputWords.insert(inputWords.end(), testWords.begin(), testWords.end());
    testWords.clear();
    int tAmount = (int) inputWords.size()/nTestShares;
    for (int i = 0; i< tAmount; i++) {
        testWords.push_back(inputWords[0]);
        inputWords.erase(inputWords.begin());
    }
    return true;
}