//
// Created by henrique on 02/11/2020.
//

#include <catch2/catch.hpp>
#include "../source/Grammar.h"
using namespace std;

vector<Symbol> terms = {Symbol("0", 0, true, false), Symbol("1", 1, true,false), Symbol("2", 2, true, false)};

vector<vector<Symbol>> words = {{Symbol("1", 1, true, false), Symbol("1", 1, true, false), Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                {Symbol("0", 0, true, false), Symbol("0", 0, true,false),  Symbol("1", 1, true,false),  Symbol("1", 1,true,false)},
                                {Symbol("0", 0, true, false), Symbol("1", 1,true,false)}};

vector<vector<Symbol>> words2 = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                 {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                 {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                 {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                 {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},{Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},{Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},{Symbol("0", 0, true, false), Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false), Symbol("1", 1, true, false)},{Symbol("0", 0, true, false), Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 0, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},{Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},{Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false)},

                                 {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},

                                 {Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("0", 0, true, false)},
                                 {Symbol("0", 0, true, false),Symbol("1", 1, true, false),Symbol("1", 1, true, false),Symbol("0", 0, true, false),Symbol("1", 1, true, false)}};


SCENARIO( "Training N-Gram", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words2, 4, make_pair(0,0));
    g.train(0,10);
    g.printGrammar();
}

SCENARIO( "Training BW", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words2, 3, make_pair(0,0));
    g.train(0, 10);
    g.printGrammar();
}

SCENARIO( "Training CGS", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words2, 3, make_pair(0,0));
    g.train(1, 10);
    g.printGrammar();
}
SCENARIO( "Training ALERGIA", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words2, 3, make_pair(0,0));
    g.train(2,10);
    g.printGrammar();
}

SCENARIO( "Training IO", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words, 2, make_pair(0,0));
    g.train(0, 10);
    g.printGrammar();
}

SCENARIO( "Training MHPCFG", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words, 2, make_pair(0,0));
    g.train(1, 10);
    g.printGrammar();
}

SCENARIO( "Training GSPCFG", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words, 2, make_pair(0,0));
    g.train(2, 10);
    g.printGrammar();
}

SCENARIO( "Training MHPCSG", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words, 1, make_pair(1,0));
    g.train(0, 10);
    g.printGrammar();
}

SCENARIO( "Training GSPCSG", "[testGrammar.cpp]" ) {
    Grammar g = Grammar(terms, 3, words, 1, make_pair(1,0));
    g.train(1, 10);
    g.printGrammar();
}