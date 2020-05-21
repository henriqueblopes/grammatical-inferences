#include <iostream>
#include "Symbol.h"
#include "Rule.h"
#include "Grammar.h"
#include <vector>
#include "peglib.h"
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <stdio.h>
#include <stdarg.h>
#include "InputWords.h"
#include <string.h>


using namespace std;


int main(int argc, char** argv) {

    InputWords iw = InputWords(false, 5);
    iw.readWords();
    iw.selectTrainingWords(172, false);
    vector<Symbol> terminals = iw.generateTerminals();
    int nInput = 10;
    vector<vector<Symbol>> iWordsLim;
    for (int i = 0; i < nInput; i++)
        iWordsLim.push_back(iw.inputWords[i]);

    pair<int,int> contextSize;

    contextSize = make_pair(0,0);
    Grammar g = Grammar(terminals, contextSize, 5, iWordsLim, 1);
    g.train(3,100);
    double p1 = g.perplexity(iw.testWords,true);
    cout<<endl <<"Final Perplexity: " <<p1 << endl;
    g.~Grammar();

    contextSize = make_pair(1,0);
    Grammar gkl = Grammar(terminals, contextSize, 5, iWordsLim, 2);
    gkl.train(4,100);
    double pkl1 = gkl.perplexity(iw.testWords,true);
    cout<< "Final Perplexity: " <<pkl1<< endl;
    gkl.~Grammar();

/*
    contextSize = make_pair(2,0);
    Grammar gkl2 = Grammar(terminals, contextSize, 10, iWordsLim, 2);
    gkl2.train(2,0);
    double pkl2 = gkl2.perplexity(iw.testWords,true);
    cout<< "Final Perplexity: " <<gkl2.perplexity(iw.testWords,true);
*/

    cout<< "Final Perplexity: g: " << p1 << " gkl: " << pkl1 << " gkl2: ";
    //Grammar g2 = Grammar(terminals, t1, contextSize, 4, 2);


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

/*std::vector<Symbol> words[] = {{Symbol("0", 0, true, false), Symbol("1", 1, true,false), Symbol("0", 0, true,false), Symbol("1", 1,true,false)},
                            {Symbol("0", 0, true, false), Symbol("0", 0, true,false), Symbol("1", 1, true,false), Symbol("1", 1,true,false)},
                            {Symbol("0", 0, true, false), Symbol("1", 1,true,false)},
                            {Symbol("0", 0, true, false), Symbol("0", 0, true,false), Symbol("1", 1, true,false), Symbol("1", 1,true,false)},
};*/