//
// Created by henrique on 02/04/19.
//

#include <iostream>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <chrono>
#include <utility>
#include "Grammar.h"

using namespace std;

Grammar::Grammar(std::vector<Symbol> terminals, pair<int,int> contextSize, int nNonTerminals,
                 std::vector<std::vector<Symbol>> words, int trainingMethod): terminals(std::move(terminals)), contextSize(std::move(contextSize)), nNonTerminals(nNonTerminals), words(std::move(words)) {
    nTerminals = terminals.size();
    generateNonTermnals();
    generateRulesCNF();
    start = rules[0].left.front();
}


Grammar::Grammar(std::vector<Symbol> terminals, std::vector<Symbol> nonterminals, std::vector<Rule> rules,
                 Symbol start, int nNonTerminals) : terminals(std::move(terminals)),
                                                    nonterminals(std::move(nonterminals)), rules(std::move(rules)),
                                                    start(std::move(start)), nNonTerminals(nNonTerminals) {}

Grammar::Grammar(std::vector<Symbol> terminals, std::vector<Symbol> nonterminals,
                 std::vector<Rule> rules, Symbol start, int maxRightHandSize, int maxProductionRules,
                 int nTerminals, int nNonTerminals) : terminals(std::move(terminals)), nonterminals(std::move(nonterminals)), rules(std::move(rules)),
                                                      start(std::move(start)), nTerminals(nTerminals),
                                                      nNonTerminals(nNonTerminals) {}

/*Grammar::Grammar(const std::vector<Symbol> &terminals, std::pair<int, int> contextSize,
                 int nNonTerminals, std::vector<std::vector<Symbol>> words, int trainingM) : terminals(terminals), contextSize(std::move(std::move(contextSize))), nNonTerminals(nNonTerminals), words(std::move(words)) {
    nTerminals = terminals.size();
    generateNonTermnals();
    generateRulesCNF();
    start = rules[0].left.front();
}*/

Grammar::Grammar(const std::vector<Symbol> &terminals, int nNonTerminals, std::vector<std::vector<Symbol>> words,
                 int type) : terminals(terminals), nNonTerminals(nNonTerminals), words(std::move(std::move(words))), type(type) {
    nTerminals = terminals.size();
    if (type == 4 ) {
        generateNGramNonTerminals();
        generateNGramRules();
    } else if (type >1){
        generateNonTermnals();
        if(type == 3)
            generateRulesRegular();
        else if (type == 2) {
            contextSize = make_pair(0,0);
            generateRulesCNF();
            start = rules[0].left.front();
        }
    } else {
        contextSize = make_pair(1,0);
        nTerminals = terminals.size();
        generateNonTermnals();
        generateRulesCNF();
        start = rules[0].left.front();
    }

}

Grammar::Grammar(const vector<Symbol> &terminals, std::vector<std::vector<Symbol>> words, int type,
                 int nNonTerminals)
        : terminals(terminals), words(std::move(words)), type(type), nNonTerminals(nNonTerminals) {
    nTerminals = terminals.size();
}


void Grammar::printGrammar() {
    std::cout << "Terminals: " << std::endl;
    std::vector<Symbol>::iterator it;
    for(it = terminals.begin(); it != terminals.end(); it++) {
        std::cout << "\t";
        std::cout <<  (*it).name << std::endl;
    }

    std::cout << std::endl;

    std::cout << "NonTerminals: " << std::endl;
    for(it = nonterminals.begin(); it != nonterminals.end(); it++) {
        std::cout << "\t";
        std::cout <<  (*it).name <<std::endl;
    }

    std::cout << std::endl;

    std::vector<Rule>::iterator itr;
    std::cout << "Rules: " << std::endl;
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        std::cout << "\t";
        (*itr).printRule();
    }

    std::cout << std::endl;
}

std::string Grammar::grammarToStr() {
    std::string g;
    g.append("Terminals: ");
    g.append("\n");
    std::vector<Symbol>::iterator it;
    for(it = terminals.begin(); it != terminals.end(); it++) {
        g.append("\t");
        g.append((*it).name);
        g.append("\n");
    }

    g.append("\n");

    g.append("NonTerminals: ");
    g.append("\n");
    for(it = nonterminals.begin(); it != nonterminals.end(); it++) {
        g.append("\t");
        g.append((*it).name);
        g.append("\n");
    }
    g.append("\n");

    std::vector<Rule>::iterator itr;
    g.append("Rules: ");
    g.append("\n");
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        g.append("\t");
        g.append((*itr).ruleToStrLALR());
    }
    g.append("\n");
    return g;
}

void Grammar::printRules() {
    std::vector<Rule>::iterator itr;
    std::cout << "Rules: " << std::endl;
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        std::cout << "\t";
        (*itr).printRule();
    }

    std::cout << std::endl;
}

std::string Grammar::rulesToString () {
    std::string rulesStr;
    std::vector<Rule>::iterator itr;
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        rulesStr += R"(
                    )" +(*itr).ruleToStr();
    }
    return rulesStr;
}

void Grammar::generateNonTermnals() {

    for (int i = 0; i< nNonTerminals; i++) {
        nonterminals.emplace_back("NT" + std::to_string(i), i, false, false);
    }
}



void Grammar::generatePermutation(std::vector<std::vector<Symbol>> & permutations, std::vector<Symbol> symbols, int size, std::vector<Symbol> word, bool context) {
    if (size == 0 )
        permutations.push_back(word);
    else {

        std::vector<Symbol>::iterator itSymbol;
        for (itSymbol = symbols.begin(); itSymbol != symbols.end(); itSymbol++) {
            Symbol aux = (*itSymbol).clone();
            aux.context = context;
            word.push_back(aux);
            generatePermutation(permutations, symbols, size-1, word, context);
            word.pop_back();
        }
    }
}

void Grammar::train(int algorithm, int iterations) {
    if (algorithm == 1) {
        contextSize.first = contextSize.second = 0;
        metropolisHastingsPCFG(iterations);
    }
    else if (algorithm == 3) {
        contextSize.first = contextSize.second = 0;
        gibbsSamplingPCFG(iterations);
    }
    else if (algorithm == 4) {
        gibbsSamplingPCSG(iterations);
    }
    else
        metropolisHastingsPCSG(iterations);
    //printGrammar();
}


double*** Grammar::CYKProb(const std::string& w) {
    //std::cout <<"IT for " << w << std::endl;
    auto ***p = new double**[nonterminals.size()];
    for (int i = 0; i < nonterminals.size(); i++) {
        p[i] = new double *[w.size()];
        for (int j = 0; j < w.size(); j++)
            p[i][j] = new double[w.size()];
    }

    for (int i = 0; i < nonterminals.size(); i++)
        for (int j = 0; j < w.size(); j++)
            for (int k = 0; k < w.size(); k++)
                p[i][j][k] = 0.0;

    std::vector<Rule>::iterator itRule;
    for (int i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol>::iterator itLeft;
            Symbol nonterminal =  Symbol("", 0,false);
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context)
                    nonterminal = (*itLeft).clone();
                break;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                        if(itRightS->name == w.substr(i,1)) {
                            p[nonterminal.id][i][i] = itRight->second.first;
                            //std::cout<< "IT["<<nonterminal.id<<"]["<<i<<"]["<<i<<"] = " << itRight->second.first << " probR" << std::endl;
                            break;
                        }
                    }
                }
            }

        }
    }
    for (int i = 1; i < w.size(); i++) {
        for (int j = 0; j < w.size()-i; j++) {
            for (int k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol>::iterator itLeft;
                    Symbol nonterminal =  Symbol("", 0,false);
                    for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                        if (!itLeft->terminal && !itLeft->context)
                            nonterminal = (*itLeft).clone();
                        break;
                    }
                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (!itRightS->terminal && !itRightS->context) {
                                double bInside = p[(*itRightS).id][j][j+k];
                                itRightS++;
                                double cInside = p[(*itRightS).id][j+k+1][j+i];
                                p[nonterminal.id][j][i+j] = p[nonterminal.id][j][i+j] + bInside*cInside*(*itRight).second.first;
                                /*std::cout<< "IT["<<nonterminal.id<<"]["<<j<<"]["<<i+j<<"] = "
                                << " IT["<<nonterminal.id<<"]["<<j<<"]["<<i+j<<"] + IT[" <<(*itRightS).id<<"]["<<j<<"]["<<j+k<<"]" <<
                                " * IT[" <<(*itRightS).id<<"]["<<j+k+1<<"]["<<j+i<<"]" << " * " <<  "ProbR = "
                                << p[nonterminal.id][j][i+j] << std::endl;*/
                                //std::cout << p[nonterminal.id][j][i+j]  << " + " << p[(*itRightS).id][j][j+k] << " * " << p[(*itRightS).id][j+k+1][j+i] << " * " << (*itRight).second.first << std::endl;
                                break;
                            }
                        }
                    }

                }
            }
        }
    }

    //printInsideTable(p,w.size());
    return p;
}

void Grammar::printInsideTable(double ***p, int wSize) const {
    cout << "IutsideTable" << endl;
    for (int i = 0; i< nonterminals.size(); i++) {
        std::cout <<"Nonterminal: " << i << "\n";
        for (int j = 0; j < wSize; j++) {
            for (int k = 0; k < wSize; k++) {
                std::cout << p[i][j][k] << " ";
            }
            std::cout << "\n";
        }
    }
}

void Grammar::printOutsideTable(const std::vector<std::vector<std::vector<double>>>& p) {
    int i = 0;
    cout << "OutsideTable" << endl;
    for (const auto& a: p) {
        cout << "Nonterminal: " << i << endl;
        i++;
        for (const auto& b: a) {
            for (auto c: b) {
                std::cout << c << " ";
            }
            std::cout << "\n";
        }
    }
    std::cout << "\n";
}

void Grammar::printInsideTableKL(double ****p, int wSize) {

    int sizeLeftContext = 0;
    int sizeRightContext = 0;
    for (int i = 0; i <= contextSize.first; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();
    }

    for (int i = 0; i <= contextSize.second; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, nonterminals, i, word, true);
        sizeRightContext += permutationsTerminals.size();
    }


    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++) {
        std::cout <<"Nonterminal: " << i << "\n";

        for (int j = 0; j < sizeRightContext; j++) {
            std::cout << "Context: " << j << "\n";
            for (int k = 0; k < wSize; k++) {
                for (int l = 0; l < wSize; l++) {
                    std::cout << p[i][l][k][j] << " ";
                }
            std::cout << "\n";
            }
        }
        std::cout << "\n";
    }
}

void Grammar::sampleParseTree(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, const std::string& w, double ***insideTable, unsigned int i, unsigned int k) {
    unsigned int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    for (int l = 0; l < jRange; l++) {
                        unsigned int indexB;
                        indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[(*itSymbol).id][i][i + l];
                        itSymbol++;
                        double pCjk = insideTable[(*itSymbol).id][i + l + 1][k];
                        unsigned int indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) / insideTable[(r.left[0].id)][i][k];
                        itSymbol--;
                        /*std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first << " * " <<
                        pBij << " ("<< "IT["<<(r.left[0].id)<<"]["<<i<<"]["<<i+l<<"])"<< " * " <<
                        pCjk << " ("<< "IT["<<(r.left[0].id)<<"]["<<i+l+1<<"]["<<k<<"])" << " / " << insideTable[(r.left[0].id)][i][k]
                        <<"IT["<<(r.left[0].id)<<"]["<<i<<"]["<<k<<"]"<<std::endl;*/
                    }
                    break;
                }

            }
        }
        int l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //double p = (double) rand()/RAND_MAX;
        //std::cout << "Free: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        unsigned int j = 0;
        delete[]jbc;
        //std::vector<Symbol> rightProduction;
        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    if ((*itSymbol).id == l / (nNonTerminals*jRange)) {
                        itSymbol++;
                        if(itSymbol == (*itRight).first.end())
                            break;
                        if((*itSymbol).id == (l%(nNonTerminals*jRange)) / jRange) {
                            itSymbol--;
                            rightProduction.first.push_back((*itSymbol));
                            j = l - nNonTerminals*jRange* (*itSymbol).id;
                            itSymbol++;
                            if (itSymbol == (*itRight).first.end())
                                std::cout << "Error\n";
                            rightProduction.first.push_back((*itSymbol));
                            if ((*itSymbol).id > 0)
                                j = j % (*itSymbol).id;
                            rightProduction.second.first = (*itRight).second.first;
                            rightProduction.second.second = (*itRight).second.second;
                            break;
                        } else {
                            itSymbol--;
                            break;
                        }
                    }
                }
            }
        }
        //TO DO:: Dar um jeito de não hardecodá o left[0]
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        sampleParseTree(vr, rules[production.second.first[0].id], w, insideTable, i, i+j);
        sampleParseTree(vr, rules[production.second.first[1].id], w, insideTable, i+j+1, k);

    } else {
        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w.substr(i,jRange+1);
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            //TO DO:: Dar um jeito de não hardecodá o left[0]
            if ((*itRight).first[0].name == s) {
                rightProduction.first.push_back((*itRight).first[0]);
                rightProduction.second.first = (*itRight).second.first;
                rightProduction.second.second = (*itRight).second.second;
            }
        }
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        return;
    }

}
void Grammar::calculateNewThetaVecOpt(const std::vector<Symbol>& w, int i) {
    pTiMinus1Frequence(w, i);
    std::vector<Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        double nonterminalTotal = 0;
        std::vector<std::pair<double, int>>::iterator iTrFrequence;
        for (iTrFrequence = (*itRule).ruleFrequence.begin(); iTrFrequence != (*itRule).ruleFrequence.end(); iTrFrequence++) {
            nonterminalTotal += iTrFrequence->second + iTrFrequence->first;
        }
        iTrFrequence = itRule->ruleFrequence.begin();
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            (*itRight).second.first = ((*iTrFrequence).second + (*itRight).second.second )/nonterminalTotal;
            iTrFrequence++;
        }
    }
    pTiMinus1PlusFrequence(w, i);
}


void Grammar::calculateNewTheta(const std::string& w) {
    std::vector<Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::vector<int> nonterminalCounts;
        double nonterminalTotal = 0;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production;
            production.first = (*itRule).left;
            production.second.first = (*itRight).first;
            int productionCount = calculateProductonCounts(production, w);
            nonterminalCounts.push_back(productionCount);
            nonterminalTotal += productionCount + (*itRight).second.second;
            //std::cout << w << " - PC: " << productionCount << " NTCount: " << nonterminalTotal << " ";
            //printOneProduction(production);
        }
        std::vector<int>::iterator itInt;
        itInt = nonterminalCounts.begin();
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            (*itRight).second.first = ((*itInt) + (*itRight).second.second )/nonterminalTotal;
            itInt++;
        }
        //(*itRule).printRule();

    }
}

int Grammar::calculateProductonCounts(const std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>& production, const std::string& w) {
    int productionCount = 0;

    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>>>>::iterator itParseTree;
    for (itParseTree = parseTrees.begin(); itParseTree != parseTrees.end(); itParseTree++) {
        if (w != (*itParseTree).first) {
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>::iterator itProduction;
            for (itProduction = (*itParseTree).second.begin(); itProduction != (*itParseTree).second.end(); itProduction++) {
                if(equalProductions((*itProduction),production)) {
                    productionCount++;
                }
            }
        }
    }
    /*if (productionCount){
        std::cout << productionCount << ": ";
        printOneProduction(production);
    }*/

    return productionCount;
}

std::vector<std::pair<double, int>> Grammar::calculateRuleFrequence(Rule r, const std::string& w) {
    std::vector<std::pair<double, int>> ruleFrequence;
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production;
        production.first = r.left;
        production.second.first = (*itRight).first;
        std::pair<double, int> productionFrequence;
        productionFrequence.first = (*itRight).second.second;
        productionFrequence.second = calculateProductonCounts(production, w);
        ruleFrequence.push_back(productionFrequence);
    }
    return ruleFrequence;
}

double Grammar::cConstant(std::vector<std::pair<double, int>> ruleFrequence) {
    double numerator = 1.0;
    double denominator = 0.0;
    std::vector<std::pair<double, int>>::iterator itFrequence;
    for (itFrequence = ruleFrequence.begin(); itFrequence != ruleFrequence.end(); itFrequence++) {
        numerator = numerator * tgamma((*itFrequence).first + (*itFrequence).second);
        denominator = denominator + (*itFrequence).first + (*itFrequence).second;
    }
    return numerator/tgamma(denominator);
}

void Grammar::pTiTiMinus1 (const std:: string& w) {
    std::vector<Rule>::iterator itRule;
    double result = 1;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        result *=  cConstant(calculateRuleFrequence((*itRule), ""))/cConstant(calculateRuleFrequence((*itRule), w));

    }
}

double Grammar::probTree(std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> tree) {
    std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>::iterator itTree;
    double prob = 1.0;
    for (itTree = tree.begin(); itTree != tree.end(); itTree++) {
        prob *= (*itTree).second.second.first;
    }
    return prob;
}

bool Grammar::equalProductions (std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> prodA, std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> prodB) {
    std::vector<Symbol>::iterator itA;
    std::vector<Symbol>::iterator itB;
    bool flag = true;
    itB = prodB.first.begin();
    for (itA = prodA.first.begin(); itA != prodA.first.end(); itA++) {
        if(!(*itA).equalSymbol((*itB))) {
            flag = false;
            break;
        }
        itB++;
        if(itB == prodB.first.end() && (itA+1) != prodA.first.end()) {
            flag = false;
            break;
        }

    }
    if(itB != prodB.first.end())
        flag = false;
    if(flag) {
        itB = prodB.second.first.begin();
        for (itA = prodA.second.first.begin(); itA != prodA.second.first.end(); itA++) {
            //retirar a condição abaixo caso dê merda em outros métodos
            if (itB == prodB.second.first.end()) {
                flag = false;
                break;
            }
            if(!(*itA).equalSymbol((*itB))) {
                flag = false;
                break;
            }
            itB++;
            if(itB == prodB.second.first.end() && (itA+1) != prodA.second.first.end()) {
                flag = false;
                break;
            }
        }
        if(itB != prodB.second.first.end())
            flag = false;
    }
    return flag;
}

void Grammar::updateParseTressTheta () {
    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>>>::iterator itParseTrees;

    for (itParseTrees = parseTrees.begin(); itParseTrees != parseTrees.end(); itParseTrees++) {
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>::iterator itProduction;
        for (itProduction = (*itParseTrees).second.begin(); itProduction != (*itParseTrees).second.end(); itProduction++) {
            std::vector<Rule>::iterator itRule;
            bool flagProduction = false;
            for(itRule = rules.begin(); itRule != rules.end(); itRule++) {
                std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                    std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production;
                    production.first = (*itRule).left;
                    production.second = (*itRight);
                    if (equalProductions((*itProduction), production)) {
                        (*itProduction) = production;
                        flagProduction = true;
                        break;
                    }
                }
                if(flagProduction)
                    break;
            }
        }

    }
}

void Grammar::printProduction(std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> tree) {
    std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>::iterator  itProduction;
    std::cout << "\nProductions:";
    for (itProduction = tree.begin(); itProduction != tree.end(); itProduction++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itProduction).first.begin(); itSymbol != (*itProduction).first.end(); itSymbol++)
                std::cout << (*itSymbol).name << " ";
            std::cout << " -> ";
            for (itSymbol = (*itProduction).second.first.begin(); itSymbol != (*itProduction).second.first.end(); itSymbol++)
                std::cout << (*itSymbol).name << " ";
            std::cout << "\n";
    }
}

void Grammar::printOneProduction(std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> prd) {
    std::vector<Symbol>::iterator itSymbol;
    for (itSymbol = prd.first.begin(); itSymbol != prd.first.end(); itSymbol++)
        std::cout << (*itSymbol).name << " ";
    std::cout << " -> ";
    for (itSymbol = (prd).second.first.begin(); itSymbol != (prd).second.first.end(); itSymbol++)
        std::cout << (*itSymbol).name << " ";
    std::cout << "\n";
}


void Grammar::metropolisHastingsPCFG(int iterations) {
    //srand((unsigned) time(nullptr));
    int sentencesLength = words.size();
    for (int i = 0; i< sentencesLength; i++) {
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
        double ***insideTable;
        insideTable = CYKProbVec(words[i]);
        sampleParseTreeVec(vProds, rules[0], words[i], insideTable, 0, words[i].size()-1);
        freeInsideTable(insideTable, words[i].size());
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTree;
        pairTree.first = words[i];
        pairTree.second = vProds;
        parseTreesVec.push_back(pairTree);
    }
    int countAcceptedTress = 0;
    bool acceptTree = true;
    for (int j = 0; j < iterations; j++) {
        double ***insideTable;

        std::mt19937 mt(rd());
        std::uniform_int_distribution<int> dist(0, sentencesLength-1);
        int p = dist(mt);
        int i = p;
        //int i = rand()%sentencesLength;

        //printInsideTable(insideTable, sentences[i].size());
        if (j%(iterations/10) == 0) {
            //std::pair<double,  double> pMpP = perplexity(words, true);
            std::cout << "iteration: " << j << " Tree " << i <<" Accepted Tress: " << countAcceptedTress << std::endl;
        }

        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
        //std::cout << "FREE" << std::endl;
        //if (acceptTree) {
            calculateNewThetaVec(words[i]);
            updateParseTressTheta();
            acceptTree = false;
        //}
        insideTable = CYKProbVec(words[i]);
        //printInsideTable(insideTable, sentences[i].size());
        vProds = parseTreesVec[i].second;
        double pTiWi = probTree(vProds)/insideTable[0][0][words[i].size()-1] ;
        double pTiTiMinis1 = pTiTiMinus1Vec(words[i]);
        vProds.clear();

        sampleParseTreeVec(vProds, rules[0], words[i], insideTable, 0, words[i].size()-1);
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> tp;
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> tpaux;
        tpaux = parseTreesVec[i];
        tp.first = parseTreesVec[i].first,
                tp.second = vProds;
        parseTreesVec[i] = tp;
        double pTiLineWi =  probTree(vProds) / insideTable[0][0][words[i].size()-1];
        double pTiLineTiMinis1 = pTiTiMinus1Vec(words[i]);
        double funcA = std::min(1.0,(pTiLineTiMinis1*pTiWi)/(pTiTiMinis1*pTiLineWi));
        double r = dist(mt);
        freeInsideTable(insideTable, words[i].size());

        if (funcA >= r)
            parseTreesVec[i] = tpaux;
        else {
            acceptTree = true;
            countAcceptedTress++;
            /*std::pair<double,  double> pMpP = perplexity(words, true);
            std::cout << "  Accepted tree in iteration: " << j << " - PerplexityM: " << pMpP.first << " = PerplexityP: " << pMpP.second<< std::endl;
            std::cout << "    Pti'ti-1 = " << pTiLineTiMinis1 << " Ptiti-1 = "  << pTiTiMinis1<< " Ptiwi = " <<  pTiWi  << " Pti'wi = " << pTiLineWi << std::endl;
            std::cout << "    ratio = " << (pTiLineTiMinis1*pTiWi)/(pTiTiMinis1*pTiLineWi) << ", funcA = " << funcA << ", r = " << r << std::endl;*/
        }
    }
    //printGrammar();
}

void Grammar::metropolisHastingsPCSG(int iterations) {
    std::chrono::duration<double> totalTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> insideTTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> insideTTimeMet = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> iterationTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> newThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> updateThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> perplexityTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> remainingTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    //srand((unsigned) time(nullptr));
    int sentencesLength = words.size();
    for (int i = 0; i< sentencesLength; i++) {
        if (i%(sentencesLength/10) ==0)
            std::cout << 100*(i/(1.0*sentencesLength))<< "% of trees parsed" << std::endl;
        auto startIt = std::chrono::system_clock::now();
        actualProduction.clear();
        actualProduction.push_back(nonterminals[0]);
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
        double ****insideTableKL ;
        auto startInsideTable = std::chrono::system_clock::now();
        insideTableKL = CYKProbKLVec(words[i]);
        insideTTime += std::chrono::system_clock::now() - startInsideTable;
        sampleParseTreeKLVec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size()-1);
        freeInsideTableKL(insideTableKL, words[i].size());
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTreeVec;
        pairTreeVec.first = words[i];
        pairTreeVec.second = vProds;
        parseTreesVec.push_back(pairTreeVec);
        iterationTime += std::chrono::system_clock::now() - startIt;
    }
    int countAcceptedTress = 0;
    iterationTime = std::chrono::system_clock::now() - std::chrono::system_clock::now();
    for (int j = 0; j < iterations; j++) {
        auto startIt = std::chrono::system_clock::now();
        std::mt19937 mt(rd());
        std::uniform_int_distribution<int> dist(0, sentencesLength-1);
        int p = dist(mt);
        int i = p;
        double ****insideTableKL;
        double ****iTableN;
        actualProduction.clear();
        actualProduction.push_back(nonterminals[0]);
        if (j%(iterations/5) == 0) {
            //auto startPerplexityTime = std::chrono::system_clock::now();
            //std::pair<double,  double> pMpP = perplexityKL(words, true);
            //perplexityTime += std::chrono::system_clock::now() - startPerplexityTime;
            //std::cout << "   iteration: " << j << " Tree " << i <<" - Perplexity: "<< pMpP.first << " - PerplexityN: "<< pMpP.second <<" Accepted Tress: " << countAcceptedTress << " PTime: "<< perplexityTime.count() <<std::endl;
            std::cout << "   iteration: " << j << " Tree " << i << " Accepted Tress: " << countAcceptedTress << std::endl;
        }
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;

        //std::cout << "SENSITIVE" << std::endl;
        auto startNewTheta = std::chrono::system_clock::now();
        calculateNewThetaVecOpt(words[i], i);
        newThetaTime += std::chrono::system_clock::now() - startNewTheta;

        updateParseTressTheta();

        auto startInsideTable = std::chrono::system_clock::now();
        insideTableKL = CYKProbKLVec(words[i]);
        insideTTimeMet += std::chrono::system_clock::now() - startInsideTable;

        auto startRemaining = std::chrono::system_clock::now();
        vProds = parseTreesVec[i].second;

        double pTiWi = probTree(vProds)/ insideTableKL[0][0][words[i].size()-1][0] ;

        auto startUpdateTheta = std::chrono::system_clock::now();
        double pTiTiMinis1 = pTiTiMinus1VecOpt(words[i], i);

        updateThetaTime += std::chrono::system_clock::now() - startUpdateTheta;
        vProds.clear();


        sampleParseTreeKLVec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size()-1);
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> tp;
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> tpaux;
        tpaux = parseTreesVec[i];
        tp.first = parseTreesVec[i].first,
        tp.second = vProds;
        parseTreesVec[i] = tp;
        double pTiLineWi =  probTree(vProds) / insideTableKL[0][0][words[i].size()-1][0] ;

        double pTiLineTiMinis1 = pTiTiMinus1VecOpt(words[i], i);
        double funcA = std::min(1.0,(pTiLineTiMinis1*pTiWi)/(pTiTiMinis1*pTiLineWi));
        freeInsideTableKL(insideTableKL, words[i].size());
        p = dist(mt);
        if (funcA >= p) {
            parseTreesVec[i] = tpaux;
            pTiMinus1PlusFrequence(words[i],i);
        }
        else{
            countAcceptedTress++;
        }
        remainingTime += std::chrono::system_clock::now() - startRemaining;
        totalTime += std::chrono::system_clock::now() - startIt;
    }
    std::cout << "   iteration: " << iterations  << " Accepted Tress: " << countAcceptedTress << std::endl;
    //printGrammar();
}

void Grammar::generateRulesCNF() {

    for (int i = 0; i <= contextSize.first; i++) {
        for (int j = 0; j <= contextSize.second; j++) {
            std::vector<std::vector<Symbol>> permutationsTerminals;
            std::vector<std::vector<Symbol>> permutationsAll;
            std::vector<std::vector<Symbol>> permutationsNonterminals;
            std::vector<Symbol> word;
            std::vector<Symbol> allSymbols;
            allSymbols.insert(allSymbols.end(), terminals.begin(), terminals.end());
            allSymbols.insert(allSymbols.end(), nonterminals.begin(), nonterminals.end());
            generatePermutation(permutationsTerminals, terminals, i, word, true);
            generatePermutation(permutationsAll, nonterminals, j, word, true);
            generatePermutation(permutationsNonterminals, nonterminals, 2, word, false);
            std::vector<std::vector<Symbol>>::iterator itPermutationsTerminals;
            std::vector<std::vector<Symbol>>::iterator itPermutationsAll;
            std::vector<std::vector<Symbol>>::iterator itPermutationsNonterminals;
            std::vector<Symbol>::iterator itNonterminals;
            std::vector<Symbol>::iterator itTerminals;
            for( itPermutationsTerminals = permutationsTerminals.begin(); itPermutationsTerminals != permutationsTerminals.end(); itPermutationsTerminals++) {
                for (itNonterminals = nonterminals.begin(); itNonterminals != nonterminals.end(); itNonterminals++) {
                    for (itPermutationsAll = permutationsAll.begin(); itPermutationsAll != permutationsAll.end(); itPermutationsAll++) {
                        std::vector<Symbol> leftHandSide;
                        leftHandSide.insert(leftHandSide.end(), (*itPermutationsTerminals).begin(), (*itPermutationsTerminals).end());
                        leftHandSide.push_back(*itNonterminals);
                        std::vector<std::pair<std::vector<Symbol>, std::pair<double, double>>> rulesByLeft;
                        leftHandSide.insert(leftHandSide.end(), (*itPermutationsAll).begin(), (*itPermutationsAll).end());
                        for (itPermutationsNonterminals = permutationsNonterminals.begin(); itPermutationsNonterminals !=  permutationsNonterminals.end(); itPermutationsNonterminals++) {
                            std::pair<std::vector<Symbol>, std::pair<double, double>> rightHandSide;
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsTerminals).begin(), (*itPermutationsTerminals).end());
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsNonterminals).begin(), (*itPermutationsNonterminals).end());
                            rightHandSide.second.first = rightHandSide.second.second = 0.0;
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsAll).begin(), (*itPermutationsAll).end());
                            rulesByLeft.push_back(rightHandSide);
                        }
                        //TO DO: Colocar o FOR abaixo dentro do anterior Unindo os pares de não terminais com terminais
                        for (itTerminals = terminals.begin(); itTerminals !=  terminals.end(); itTerminals++) {
                            std::pair<std::vector<Symbol>,std::pair<double, double>> rightHandSide;
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsTerminals).begin(), (*itPermutationsTerminals).end());
                            Symbol aux = (*itTerminals).clone();
                            aux.context = false;
                            rightHandSide.first.push_back(aux);
                            rightHandSide.second.first = rightHandSide.second.second = 0.0;
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsAll).begin(), (*itPermutationsAll).end());
                            rulesByLeft.push_back(rightHandSide);
                        }
                        Rule r = Rule(leftHandSide, rulesByLeft);
                        r.generatePiorDirichlet(ALFA);
                        r.updateProbDirichletTheta();
                        r.leftContext.insert(r.leftContext.end(), (*itPermutationsTerminals).begin(), (*itPermutationsTerminals).end());
                        r.rightContext.insert(r.rightContext.end(),(*itPermutationsAll).begin(), (*itPermutationsAll).end());
                        r.index1stNonContext = r.leftContext.size();
                        rules.push_back(r);
                    }
                }
            }
        }
    }
    std::vector<Rule>::iterator itRule;
    int j = 0;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        (*itRule).index = j;
        for (int i =0; i < (*itRule).right.size(); i++)
            (*itRule).ruleFrequence.emplace_back(ALFA, 0);
        j++;
    }


}

double ****Grammar::CYKProbKL(const std::string& w) {
    //std::cout <<"Sensitive IT for " << w << std::endl;
    int sizeLeftContext = 0;
    int sizeRightContext = 0;
    for (int i = 0; i <= contextSize.first; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();

    }
    contextAmount.first = sizeLeftContext;
    for (int i = 0; i <= contextSize.second; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, nonterminals, i, word, true);
//        rightContexts.insert(rightContexts.end(), permutationsTerminals.begin(), permutationsTerminals.end());
        sizeRightContext += permutationsTerminals.size();
    }
    contextAmount.second = sizeRightContext;
    auto ****p = new double***[nonterminals.size()*sizeLeftContext];
    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++) {
        p[i] = new double **[w.size()];
        for (int j = 0; j < w.size(); j++) {
            p[i][j] = new double *[w.size()];
            for (int k = 0; k < w.size(); k++)
                p[i][j][k] = new double[sizeRightContext];
        }
    }

    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++)
        for (int j = 0; j < w.size(); j++)
            for (int k = 0; k < w.size(); k++)
                for (int l = 0; l < sizeRightContext; l++)
                    p[i][j][k][l] = 0.0;

    std::vector<Rule>::iterator itRule;
    for (int i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol>::iterator itLeft;
            Symbol nonterminal =  Symbol("", 0,false);
            std::vector<Symbol> leftContext;
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context) {
                    nonterminal = (*itLeft).clone();
                    itLeft++;
                    break;
                }
                else {
                    leftContext.push_back((*itLeft));
                }
            }
            std::vector<Symbol> rightContext;
            while (itLeft != (*itRule).left.end() ) {
                rightContext.push_back((*itLeft));
                itLeft++;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                        if(itRightS->name == w.substr(i,1)) {
                            p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][i][i][convertContextToID(1,rightContext)] = itRight->second.first;
                            //std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<i<<"]["<<i<<"]["<< convertContextToID(1,rightContext)<<"] = ProbR"
                                     /* <<(*itRight).second.first << std::endl;*/
                            break;
                        }
                    }
                }
            }

        }
    }
    for (int i = 1; i < w.size(); i++) {
        for (int j = 0; j < w.size()-i; j++) {
            for (int k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol>::iterator itLeft;
                    Symbol nonterminal =  Symbol("", 0,false);
                    std::vector<Symbol> leftContext;
                    for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                        if (!itLeft->terminal && !itLeft->context) {
                            nonterminal = (*itLeft).clone();
                            itLeft++;
                            break;
                        }
                        else {
                            leftContext.push_back((*itLeft));
                        }
                    }
                    std::vector<Symbol> rightContext;
                    while (itLeft != (*itRule).left.end() ) {
                        rightContext.push_back((*itLeft));
                        itLeft++;
                    }

                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (!itRightS->terminal && !itRightS->context) {
                                itRightS++;
                                std::vector<Symbol> rightContextLine;
                                if ( contextSize.second > 0) {
                                    rightContextLine.push_back((*itRightS));
                                    rightContextLine.front().context = true;
                                    if (!rightContext.empty()) {
                                        rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                                        rightContext.pop_back();
                                    }
                                }
                                itRightS--;
                                double bInside = p[sizeLeftContext*(*itRightS).id + convertContextToID(0,leftContext)][j][j+k][convertContextToID(1,rightContextLine)];
                                itRightS++;
                                double cInside = p[sizeLeftContext*(*itRightS).id + convertContextToID(0,leftContext)][j+k+1][j+i][convertContextToID(1,rightContext)];
                                p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][j][i+j][convertContextToID(1,rightContext)] =
                                        p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][j][i+j][convertContextToID(1,rightContext)] + bInside*cInside*(*itRight).second.first;

                                /*std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,rightContext)<<"] = "
                                         << "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,rightContext)<<
                                         "] + IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<j+k<<"]["<< convertContextToID(1,rightContext)<<
                                         "] * IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,leftContext)<<"]["<<j+k+1<<"]["<<j+i<<"]["<< convertContextToID(1,rightContext) << "] * " <<  "ProbR = "
                                         <<p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][j][i+j][convertContextToID(1,rightContext)]<< std::endl;*/
                                         //std::cout << p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][j][i+j][convertContextToID(1,rightContext)]  << " + " << bInside << " * " << cInside << " * " << (*itRight).second.first << std::endl;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return p;
}

int Grammar::convertContextToID(int side, std::vector<Symbol> context) const {
    int nSymbols;
    int id = 0;
    if (context.empty())
        return id;
    if (side == 0)
        nSymbols = nTerminals;
    else
        nSymbols = nNonTerminals;
    std::vector<Symbol>::iterator itSymbol;
    for (int i = 0; i < context.size(); i++)
        id += pow(nSymbols, i);
    int position = 0;
    for (itSymbol = context.begin(); itSymbol != context.end(); itSymbol++) {
        id += static_cast<int>((*itSymbol).id * pow(nSymbols, position));
        position++;
    }
    return id;
}

void Grammar::sampleParseTreeKL(
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> &vr,
        Rule r, const std::string& w, double ****insideTable, unsigned int i, unsigned int k) {
    unsigned int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    Symbol nonterminal =  Symbol("", 0,false);
    std::vector<Symbol>::iterator itLeft;
    std::vector<Symbol> leftContext;
    for (itLeft = r.left.begin(); itLeft != r.left.end(); itLeft++) {
        if (!itLeft->terminal && !itLeft->context) {
            nonterminal = (*itLeft).clone();
            itLeft++;
            break;
        }
        else {
            leftContext.push_back((*itLeft));
        }
    }
    std::vector<Symbol> rightContext;
    while (itLeft != r.left.end() ) {
        rightContext.push_back((*itLeft));
        itLeft++;
    }
    double sumJBC = 0.0;
    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    std::vector<Symbol> rightContextLine;
                    itSymbol++;
                    if ( contextSize.second > 0) {
                        rightContextLine.push_back((*itSymbol));
                        rightContextLine.front().context = true;
                        if (!rightContext.empty()) {
                            rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                            rightContext.pop_back();
                        }
                    }
                    itSymbol--;
                    for (int l = 0; l < jRange; l++) {
                        unsigned int indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                                [i][i + l][convertContextToID(1,rightContextLine)];
                        itSymbol++;
                        double pCjk = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                                [i + l + 1][k][convertContextToID(1,rightContext)];
                        unsigned int indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) /insideTable
                                [contextAmount.first*nonterminal.id+ convertContextToID(0,leftContext)]
                                [i][k][convertContextToID(1,rightContext)];
                        itSymbol--;
                        /*std::cout << "SumJBC = " <<sumJBC << " ";

                        std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first <<
                        " * " << pBij << " ("<< "IT["<<
                                      leftContext.size()*(*itSymbol).id+ convertContextToID(0,leftContext)<<"]["<<i<<"]["<<i+l<<"]["<<
                                      convertContextToID(1,rightContext)<<"])"
                        << " * " << pCjk
                                << " ("<< "IT["<<
                                leftContext.size()*(*itSymbol).id+ convertContextToID(0,leftContext)<<"]["<<i+l+1<<"]["<<k<<"]["<<
                                convertContextToID(1,rightContext)<<"])"
                        << " / " << insideTable
                        [leftContext.size()*nonterminal.id+ convertContextToID(0,leftContext)]
                        [i][k][convertContextToID(1,rightContext)] << " IT["<<
                        leftContext.size()*nonterminal.id+ convertContextToID(0,leftContext)<<"]["<<i<<"]["<<k<<"]["<<
                        convertContextToID(1,rightContext)<<"]"<<std::endl;*/


                    }
                    break;
                }
            }
        }


        int l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //std::cout << "Sensitive: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        unsigned int j = 0;

        delete[]jbc;
        Symbol rightB = Symbol("", 0,false);
        Symbol rightC = Symbol("", 0,false);

        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    if ((*itSymbol).id == l / (nNonTerminals*jRange)) {
                        itSymbol++;
                        if(itSymbol == (*itRight).first.end())
                            break;
                        if((*itSymbol).id == (l%(nNonTerminals*jRange)) / jRange) {
                            itSymbol--;
                            rightProduction.first.push_back((*itSymbol));
                            rightB = (*itSymbol).clone();
                            j = l - nNonTerminals*jRange* (*itSymbol).id;
                            itSymbol++;
                            if (itSymbol == (*itRight).first.end())
                                std::cout << "Error\n";
                            rightProduction.first.push_back((*itSymbol));
                            rightC = (*itSymbol).clone();
                            if ((*itSymbol).id > 0)
                                j = j % (*itSymbol).id;
                            rightProduction.second.first = (*itRight).second.first;
                            rightProduction.second.second = (*itRight).second.second;
                            break;
                        } else {
                            itSymbol--;
                            break;
                        }
                    }
                }
            }
        }

        rightProduction.first.insert(rightProduction.first.end(), rightContext.begin(), rightContext.end());
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        applyProduction(production);
        std::vector<Symbol> leftContext2;
        std::vector<Symbol> rightContext2;
        getActualContext(leftContext2,rightContext2 );

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();

        unsigned int lc , rc;
        std::uniform_int_distribution<int> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<int> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = ((int) rand()%(leftContext2.size()+1));
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::vector<Symbol> lhs;
        lhs.insert(lhs.end(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightB);
        lhs.insert(lhs.end(), rightContext2.begin(), rightContext2.end());
        Rule nextRule = findRuleByLHS(lhs);
        sampleParseTreeKL(vr, nextRule, w, insideTable, i, i+j);

        getActualContext(leftContext2,rightContext2 );

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();


        lc = distIntL(mt);
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();

        lhs.clear();
        lhs.insert(lhs.begin(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightC);
        lhs.insert(lhs.begin(), rightContext2.begin(), rightContext2.end());
        sampleParseTreeKL(vr, findRuleByLHS(lhs), w, insideTable, i+j+1, k);

    } else {
        std::vector<Symbol> leftContext2;
        std::vector<Symbol> rightContext2;
        getActualContext(leftContext2, rightContext2);

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();
        unsigned int lc , rc;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<int> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<int> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();

        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w.substr(i,jRange+1);
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            std::vector<Symbol>::iterator itSymbol;
            bool flagContext = false;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).name == s && !(*itSymbol).context ) {
                    rightProduction.first.push_back((*itSymbol));
                    rightProduction.second.first = (*itRight).second.first;
                    rightProduction.second.second = (*itRight).second.second;
                    flagContext = true;
                } else if (flagContext)
                    rightProduction.first.push_back((*itSymbol));
            }

        }
        //rightProduction.first.insert(rightProduction.first.end(), rightContext.begin(), rightContext.end());
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        applyProduction(production);
        return;
    }
}

void Grammar::applyProduction(
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prod) {
    std::vector<Symbol>::iterator itLeft;
    std::vector<Symbol> newProd;
    bool appliedProd = false;
    for (itLeft = prod.first.begin(); itLeft != prod.first.end(); itLeft++) {
        if (!itLeft->terminal) {
            std::vector<Symbol>::iterator itActualProd;

            for (itActualProd = actualProduction.begin(); itActualProd != actualProduction.end(); itActualProd++) {
                if ((*itActualProd).equalSymbol(*itLeft) && !appliedProd) {
                    std::vector<Symbol>::iterator  itRight;
                    for (itRight = prod.second.first.begin(); itRight != prod.second.first.end(); itRight++) {
                        if (!itRight->terminal && !itRight->context) {
                            newProd.push_back(*itRight);
                            itRight++;
                            if (itRight->context)
                                std::cout << "Maybe Error" << std::endl;
                            newProd.push_back(*(itRight));
                            appliedProd = true;
                            break;
                        }
                        if (itRight->terminal && !itRight->context) {
                            newProd.push_back(*itRight);
                            appliedProd = true;
                            break;
                        }
                    }
                } else {
                    newProd.push_back(*itActualProd);
                }
            }
        }
        if (appliedProd)
            break;
    }
    actualProduction = newProd;
}

Rule Grammar::findRuleByLHS(std::vector<Symbol> lhs) {
    std::vector<Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<Symbol>::iterator itLHS;
        auto itRuleLeft = (*itRule).left.begin();
        bool flag = true;
        for (itLHS = lhs.begin(); itLHS != lhs.end(); itLHS++) {

            if ((*itRuleLeft).name != (*itLHS).name)
                flag = false;
            if (itRuleLeft+1 != (*itRule).left.end())
                itRuleLeft++;
        }
        if (flag)
            return (*itRule);
    }
    return rules[0];
}

void Grammar::getActualContext(std::vector<Symbol> &leftContext, std::vector<Symbol> &rightContext) {
    leftContext.clear();
    rightContext.clear();
    std::vector<Symbol>::iterator itLeft;
    Symbol aux = Symbol("",0,true);
    for (itLeft = actualProduction.begin(); itLeft != actualProduction.end(); itLeft++) {
        if (!itLeft->terminal && !itLeft->context) {
            itLeft++;
            break;
        }
        else {
            aux = *itLeft;
            aux.context = true;
            aux.terminal = true;
            leftContext.push_back(aux);
        }
    }
    while (itLeft != actualProduction.end() ) {
        aux = *itLeft;
        aux.context = true;
        aux.terminal = false;
        rightContext.push_back(aux);
        itLeft++;
    }
}

void Grammar::sampleParseTreeVec(
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> & vr,
        Rule r, std::vector<Symbol> w, double ***insideTable, unsigned int i, unsigned int k) {
    unsigned int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    for (int l = 0; l < jRange; l++) {
                        unsigned int indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[(*itSymbol).id][i][i + l];
                        itSymbol++;
                        double pCjk = insideTable[(*itSymbol).id][i + l + 1][k];
                        unsigned int indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) / insideTable[(r.left[0].id)][i][k];
                        itSymbol--;
                        /*std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first << " * " <<
                        pBij << " ("<< "IT["<<(r.left[0].id)<<"]["<<i<<"]["<<i+l<<"])"<< " * " <<
                        pCjk << " ("<< "IT["<<(r.left[0].id)<<"]["<<i+l+1<<"]["<<k<<"])" << " / " << insideTable[(r.left[0].id)][i][k]
                        <<"IT["<<(r.left[0].id)<<"]["<<i<<"]["<<k<<"]"<<std::endl;*/
                    }
                    break;
                }

            }
        }
        int l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //std::cout << "Free: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        unsigned int j = 0;
        delete[]jbc;
        //std::vector<Symbol> rightProduction;
        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    if ((*itSymbol).id == l / (nNonTerminals*jRange)) {
                        itSymbol++;
                        if(itSymbol == (*itRight).first.end())
                            break;
                        if((*itSymbol).id == (l%(nNonTerminals*jRange)) / jRange) {
                            itSymbol--;
                            rightProduction.first.push_back((*itSymbol));
                            j = l - nNonTerminals*jRange* (*itSymbol).id;
                            itSymbol++;
                            if (itSymbol == (*itRight).first.end())
                                std::cout << "Error\n";
                            rightProduction.first.push_back((*itSymbol));
                            if ((*itSymbol).id > 0)
                                j = j % (*itSymbol).id;
                            rightProduction.second.first = (*itRight).second.first;
                            rightProduction.second.second = (*itRight).second.second;
                            break;
                        } else {
                            itSymbol--;
                            break;
                        }
                    }
                }
            }
        }
        //TO DO:: Dar um jeito de não hardecodá o left[0]
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        sampleParseTreeVec(vr, rules[production.second.first[0].id], w, insideTable, i, i+j);
        sampleParseTreeVec(vr, rules[production.second.first[1].id], w, insideTable, i+j+1, k);

    } else {
        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w[i].name;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            //TO DO:: Dar um jeito de não hardecodá o left[0]
            if ((*itRight).first[0].name == s) {
                rightProduction.first.push_back((*itRight).first[0]);
                rightProduction.second.first = (*itRight).second.first;
                rightProduction.second.second = (*itRight).second.second;
            }
        }
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        return;
    }

}

double ***Grammar::CYKProbVec(std::vector<Symbol> w) {
    //std::cout <<"IT for " << w << std::endl;
    auto ***p = new double**[nonterminals.size()];
    for (int i = 0; i < nonterminals.size(); i++) {
        p[i] = new double *[w.size()];
        for (int j = 0; j < w.size(); j++)
            p[i][j] = new double[w.size()];
    }

    for (int i = 0; i < nonterminals.size(); i++)
        for (int j = 0; j < w.size(); j++)
            for (int k = 0; k < w.size(); k++)
                p[i][j][k] = 0.0;

    std::vector<Rule>::iterator itRule;
    for (int i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol>::iterator itLeft;
            Symbol nonterminal =  Symbol("", 0,false);
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context)
                    nonterminal = (*itLeft).clone();
                break;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                        if(itRightS->name == w[i].name) {
                            p[nonterminal.id][i][i] = itRight->second.first;
                            //std::cout<< "IT["<<nonterminal.id<<"]["<<i<<"]["<<i<<"] = " << itRight->second.first << " probR" << std::endl;
                            break;
                        }
                    }
                }
            }

        }
    }
    for (int i = 1; i < w.size(); i++) {
        for (int j = 0; j < w.size()-i; j++) {
            for (int k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol>::iterator itLeft;
                    Symbol nonterminal =  Symbol("", 0,false);
                    for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                        if (!itLeft->terminal && !itLeft->context)
                            nonterminal = (*itLeft).clone();
                        break;
                    }
                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (!itRightS->terminal && !itRightS->context) {
                                double bInside = p[(*itRightS).id][j][j+k];
                                itRightS++;
                                double cInside = p[(*itRightS).id][j+k+1][j+i];
                                p[nonterminal.id][j][i+j] = p[nonterminal.id][j][i+j] + bInside*cInside*(*itRight).second.first;
                                /*std::cout<< "IT["<<nonterminal.id<<"]["<<j<<"]["<<i+j<<"] = "
                                << " IT["<<nonterminal.id<<"]["<<j<<"]["<<i+j<<"] + IT[" <<(*itRightS).id<<"]["<<j<<"]["<<j+k<<"]" <<
                                " * IT[" <<(*itRightS).id<<"]["<<j+k+1<<"]["<<j+i<<"]" << " * " <<  "ProbR = "
                                << p[nonterminal.id][j][i+j] << std::endl;*/
                                //std::cout << p[nonterminal.id][j][i+j]  << " + " << p[(*itRightS).id][j][j+k] << " * " << p[(*itRightS).id][j+k+1][j+i] << " * " << (*itRight).second.first << std::endl;
                                break;
                            }
                        }
                    }

                }
            }
        }
    }

    //printInsideTable(p,w.size());
    return p;
}

vector<vector<vector<double>>> Grammar::outsideTable(const std::vector<Symbol>& w, double *** insideTable) {
    vector<vector<vector<double>>>  outsideTable;
    for (const auto& nt: nonterminals) {
        outsideTable.emplace_back();
        for (int i = 0; i < w.size(); i++) {
            outsideTable[nt.id].push_back(vector<double>());
            for (int j = 0; j < w.size(); j++) {
                outsideTable[nt.id][i].push_back(0.0);
            }
        }
    }
    outsideTable[0][0][w.size()-1] = 1;
    //printOutsideTable(outsideTable);

    for (unsigned int i = w.size()-1; i >=0; i--) {
        for (int j = 0; j+i< w.size(); j++) {
            for (const auto& nt: nonterminals) {
                //cout << "B_"<<nt.id<<"_"<<j<<"_"<<j+i<<endl;

                for (auto r: rules) {
                    for (auto right: r.right) {
                        if (right.first[0].terminal)
                            break;
                        for (int k=0; k < j; k++) {

                            if (right.first[1].id ==nt.id) {
                                outsideTable[nt.id][j][j+i] += outsideTable[r.left[0].id][k][j+i]*insideTable[ right.first[0].id][k][j-1]* right.second.first;
                                /*cout << "B_" << r.left[0].id << "_" << k << "_" << j + i << " * " << "I_" << right.first[0].id << "_" << k << "_" << j - 1 << " = "
                                     <<outsideTable[r.left[0].id][k][j+i] << " * " << insideTable[ right.first[0].id][k][j-1] << " * " << right.second.first
                                     << " = " << outsideTable[nt.id][j][j+i] << endl;*/
                            }                    }
                        for (unsigned int k=w.size()-1; k > j+i; k--) {
                            if (right.first[0].id ==nt.id) {
                                outsideTable[nt.id][j][j+i] += outsideTable[r.left[0].id][j][k]*insideTable[ right.first[1].id][j+i+1][k]* right.second.first;
                                /*cout << "B_" <<r.left[0].id<<"_"<<j<<"_"<<k<< " * " << "I_"<<right.first[1].id<<"_"<<j+i+1<<"_"<<k << " = "
                                     << outsideTable[r.left[0].id][j][k] << " * " <<insideTable[ right.first[1].id][j+i+1][k] << " * " << right.second.first
                                     << " = " << outsideTable[nt.id][j][j+i] <<endl;*/
                            }
                        }
                    }
                }
                //cout << endl;
            }

        }
    }



    /*
    cout << endl << endl << "ordered outside:" <<endl;

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j<=n; j++) {
            cout << "B_"<<nt<<"_"<<i<<"_"<<j<<endl;
            for (auto r: rules) {
                for (auto right: r.right) {
                    if (right.first[0].terminal)
                        break;
                    for (int k=1; k < i; k++) {
                        if (right.first[1].id ==nt)
                            cout << "B_" <<r.left[0].id<<"_"<<k<<"_"<<j << " * " << "I_"<<right.first[0].id<<"_"<<k<<"_"<<i-1 <<endl;
                    }
                    for (int k=n; k > j; k--) {
                        if (right.first[0].id ==nt)
                            cout << "B_" <<r.left[0].id<<"_"<<i<<"_"<<k<< " * " << "I_"<<right.first[1].id<<"_"<<j+1<<"_"<<k <<endl;

                    }
                }
            }
            cout << endl;
        }
    }*/

    /*for (int i = w.size()-1; i>0 ; i++) {
        for (int j = i-1; j >0; j++) {

        }
    }*/

/*
    for (int e = 0; e<= w.size(); e++) {
        for (int j = 1; j <= w.size()-e; j++) {
            int k = w.size()+1-j+e;
            for (auto r: rules) {
                for (int l = 0; l < nNonTerminals*nNonTerminals; l++) {
                    for (int s = j; s <=k; s++) {
                        outsideTable[r.right[l].first[0].id][j][s] += r.right[l].second.first *
                                outsideTable[r.left[0].id][j][k] * insideTable[r.right[l].first[1].id][s][k-1];

                        cout << "O["<<r.right[l].first[0].id<<"]["<<j<<"]["<<s<<"] += O[" << r.left[0].id<<"]["<<j<<"]["<<k<<"] * I["
                        <<r.right[l].first[1].id<<"]["<<s<<"]["<<k-1<<"]\t\t" ;



                        outsideTable[r.right[l].first[1].id][s][k] += r.right[l].second.first *
                                outsideTable[r.left[0].id][j][k] * insideTable[r.right[l].first[0].id][j+1][s];

                        cout << "O["<<r.right[l].first[1].id<<"]["<<s<<"]["<<k<<"] += O[" << r.left[0].id<<"]["<<j<<"]["<<k<<"] * I["
                             <<r.right[l].first[0].id<<"]["<<j+1<<"]["<<s<<"]" << endl;
                    }
                }
            }
        }
    }*/
    /*for (int p = 0; p < w.size(); p++) {
        for (int q=0; q < w.size(); q++) {
            for (auto r: rules) {
                for (int l = 0; l < nNonTerminals*nNonTerminals; l++) {
                    for (int k = q+1; k <w.size();k++) {
                        outsideTable[r.right[l].first[0].id][p][q] += outsideTable[r.left[0].id][p][k] *
                                insideTable[r.right[l].first[1].id][q+1][k]*r.right[l].second.first;
                    }
                    for (int k = 0; k <p-1;k++) {

                    }
                }
            }
        }
    }*/
    /*for (int i = w.size()-1; i >0 ; i--) {
        for (int j = 0; j < w.size() - i; j++) {
            for (int k = 0; k < i; k++) {
                for (auto r: rules) {
                    for (int l = 0; l < nNonTerminals*nNonTerminals; l++) {
                        cout << "i:"<<i<<" j:"<<j<<" k:"<<k<<" - ";
                        outsideTable[r.right[l].first[0].id][j][j+k] +=
                                outsideTable[r.left[0].id][j][j+i] * insideTable[r.right[l].first[1].id][j+k][j+i] * r.right[l].second.first;

                        cout << "OT["<<r.right[l].first[0].id<<"]["<<j<<"]["<<j+k<<"] += OT["<<r.left[0].id<<"]["<<j<<"]["<<j+i<<"] * IT["
                             << r.right[l].first[1].id<<"]["<<j+k<<"]["<<j+i<<"] * p"<<r.left[0].id<<"->"<<r.right[l].first[0].id<<r.right[l].first[1].id << " = " << outsideTable[r.left[0].id][j][j+i] << "*"
                             << insideTable[r.right[l].first[1].id][j+k][j+i] <<"*" << r.right[l].second.first << " = " << outsideTable[r.right[l].first[0].id][j][j+k];

                        outsideTable[r.right[l].first[1].id][j+k][j+i] +=
                                outsideTable[r.left[0].id][j][j+i] * insideTable[r.right[l].first[0].id][j][j+k] * r.right[l].second.first;

                        cout << "\tOT["<<r.right[l].first[1].id<<"]["<<j+k<<"]["<<j+i<<"] += OT["<<r.left[0].id<<"]["<<j<<"]["<<j+i<<"] * IT["
                             <<r.right[l].first[0].id<<"]["<<j<<"]["<<j+k<<"] * p"<<r.left[0].id<<"->"<<r.right[l].first[0].id<<r.right[l].first[1].id << " = " << outsideTable[r.left[0].id][j][j+i] << "*"
                             << insideTable[r.right[l].first[0].id][j][j+k] <<"*" << r.right[l].second.first << " = " << outsideTable[r.right[l].first[0].id][j+k][j+i] << endl;
                    }
                }
            }
        }
    }*/
    return outsideTable;
}

void Grammar::calculateNewThetaVec(const std::vector<Symbol>& w) {
    std::vector<Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::vector<int> nonterminalCounts;
        double nonterminalTotal = 0;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production;
            production.first = (*itRule).left;
            production.second.first = (*itRight).first;
            int productionCount = calculateProductonCountsVec(production, w);
            nonterminalCounts.push_back(productionCount);
            nonterminalTotal += productionCount + (*itRight).second.second;
            //std::cout << w << " - PC: " << productionCount << " NTCount: " << nonterminalTotal << " ";
            //printOneProduction(production);
        }
        std::vector<int>::iterator itInt;
        itInt = nonterminalCounts.begin();
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            (*itRight).second.first = ((*itInt) + (*itRight).second.second )/nonterminalTotal;
            itInt++;
        }
        //(*itRule).printRule();

    }
}

int Grammar::calculateProductonCountsVec(
        const std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>& production,
        const std::vector<Symbol>& w) {
    int productionCount = 0;

    std::vector<std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>>>>::iterator itParseTreeVec;
    for (itParseTreeVec = parseTreesVec.begin(); itParseTreeVec != parseTreesVec.end(); itParseTreeVec++) {
        if (!equalWord(w,(*itParseTreeVec).first)) {
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>::iterator itProduction;
            for (itProduction = (*itParseTreeVec).second.begin(); itProduction != (*itParseTreeVec).second.end(); itProduction++) {
                if(equalProductions((*itProduction),production)) {
                    productionCount++;
                }
            }
        }
    }
    /*if (productionCount){
        std::cout << productionCount << ": ";
        printOneProduction(production);
    }*/

    return productionCount;
}

bool Grammar::equalWord(std::vector<Symbol> w1, std::vector<Symbol> w2) {
    std::vector<Symbol>::iterator itW1;
    auto itW2 = w2.begin();
    for (itW1 = w1.begin(); itW1 != w1.end(); itW1++) {
        if(itW2 == w2.end())
            return false;
        if (!(*itW1).equalSymbol((*itW2)))
            return false;
        itW2++;
    }
    if(itW2!= w2.end())
        return false;
    else
        return true;


}

double Grammar::pTiTiMinus1Vec(const std::vector<Symbol>& w) {
    std::vector<Rule>::iterator itRule;
    double result = 1.0;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<Symbol> a;
        result *=  cConstant(calculateRuleFrequenceVec((*itRule),  a))/cConstant(calculateRuleFrequenceVec((*itRule), w));
        //std::cout << "   PTia = " << cConstant(calculateRuleFrequenceVec((*itRule),  a)) << " Pt-ia = " << cConstant(calculateRuleFrequenceVec((*itRule), w)) << std::endl;
        //Talvez essa passagem de parametro esteja errada

    }
    return result;
}

double Grammar::pTiTiMinus1VecOpt(const std::vector<Symbol>& w, int i) {
    std::vector<Rule>::iterator itRule;
    std::vector<std::vector<std::pair<double, int>>> vecPairs;
    vecPairs.reserve(rules.size());
    for (itRule = rules.begin();itRule < rules.end(); itRule++) {
        vecPairs.push_back(itRule->ruleFrequence);
    }
    pTiMinus1Frequence(w,i);
    double result = 1.0;
    auto vecPairsIt = vecPairs.begin();
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<Symbol> a;
        result *=  cConstant(*vecPairsIt)/cConstant(itRule->ruleFrequence);
        vecPairsIt++;
    }
    return result;
}

void Grammar::pTiMinus1Frequence(const std::vector<Symbol>& w, int i) {
    if (i == -1)
        return;
    std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>::iterator pTreeIt;
    for (pTreeIt = parseTreesVec[i].second.begin(); pTreeIt != parseTreesVec[i].second.end(); pTreeIt++) {
        Rule r = findRuleByLHS(pTreeIt->first);
        if (!(*pTreeIt).second.first[r.index1stNonContext].terminal)
            rules[r.index].ruleFrequence[(*pTreeIt).second.first[r.index1stNonContext].id*nNonTerminals +
                (*pTreeIt).second.first[r.index1stNonContext+1].id].second--;
        else
            rules[r.index].ruleFrequence[nNonTerminals*nNonTerminals+ (*pTreeIt).second.first[r.index1stNonContext].id].second--;
    }
}

void Grammar::pTiMinus1PlusFrequence(const std::vector<Symbol>& w, int i) {
    if (i == -1)
        return;
    std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>::iterator pTreeIt;
    for (pTreeIt = parseTreesVec[i].second.begin(); pTreeIt != parseTreesVec[i].second.end(); pTreeIt++) {
        Rule r = findRuleByLHS(pTreeIt->first);
        if (!(*pTreeIt).second.first[r.index1stNonContext].terminal)
            rules[r.index].ruleFrequence[(*pTreeIt).second.first[r.index1stNonContext].id*nNonTerminals +
                                         (*pTreeIt).second.first[r.index1stNonContext+1].id].second++;
        else
            rules[r.index].ruleFrequence[nNonTerminals*nNonTerminals+ (*pTreeIt).second.first[r.index1stNonContext].id].second++;
    }
}

std::vector<std::pair<double, int>> Grammar::calculateRuleFrequenceVec(Rule &r, const std::vector<Symbol>& w) {
    std::vector<std::pair<double, int>> ruleFrequence;
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production;
        production.first = r.left;
        production.second.first = (*itRight).first;
        std::pair<double, int> productionFrequence;
        productionFrequence.first = (*itRight).second.second;
        productionFrequence.second = calculateProductonCountsVec(production, w);
        ruleFrequence.push_back(productionFrequence);
    }
    return ruleFrequence;
}



double ****Grammar::CYKProbKLVec(std::vector<Symbol> w) {
    int sizeLeftContext = 0;
    int sizeRightContext = 0;
    for (int i = 0; i <= contextSize.first; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();

    }
    contextAmount.first = sizeLeftContext;
    for (int i = 0; i <= contextSize.second; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, nonterminals, i, word, true);
//        rightContexts.insert(rightContexts.end(), permutationsTerminals.begin(), permutationsTerminals.end());
        sizeRightContext += permutationsTerminals.size();
    }
    contextAmount.second = sizeRightContext;
    auto ****p = new double***[nonterminals.size()*sizeLeftContext];
    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++) {
        p[i] = new double **[w.size()];
        for (int j = 0; j < w.size(); j++) {
            p[i][j] = new double *[w.size()];
            for (int k = 0; k < w.size(); k++)
                p[i][j][k] = new double[sizeRightContext];
        }
    }

    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++)
        for (int j = 0; j < w.size(); j++)
            for (int k = 0; k < w.size(); k++)
                for (int l = 0; l < sizeRightContext; l++)
                    p[i][j][k][l] = 0.0;

    std::vector<Rule>::iterator itRule;
    for (int i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol>::iterator itLeft;
            Symbol nonterminal =  Symbol("", 0,false);
            std::vector<Symbol> leftContext;
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context) {
                    nonterminal = (*itLeft).clone();
                    itLeft++;
                    break;
                }
                else {
                    leftContext.push_back((*itLeft));
                }
            }
            std::vector<Symbol> rightContext;
            while (itLeft != (*itRule).left.end() ) {
                rightContext.push_back((*itLeft));
                itLeft++;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                        if(itRightS->name == w[i].name) {
                            p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][i][i][convertContextToID(1,rightContext)] = itRight->second.first;
                            //std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<i<<"]["<<i<<"]["<< convertContextToID(1,rightContext)<<"] = ProbR"
                            /* <<(*itRight).second.first << std::endl;*/
                            break;
                        }
                    }
                }
            }

        }
    }
    for (int i = 1; i < w.size(); i++) {
        for (int j = 0; j < w.size()-i; j++) {
            for (int k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol>::iterator itLeft;
                    Symbol nonterminal =  Symbol("", 0,false);
                    std::vector<Symbol> leftContext;
                    for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                        if (!itLeft->terminal && !itLeft->context) {
                            nonterminal = (*itLeft).clone();
                            itLeft++;
                            break;
                        }
                        else {
                            leftContext.push_back((*itLeft));
                        }
                    }
                    std::vector<Symbol> rightContext;
                    while (itLeft != (*itRule).left.end() ) {
                        rightContext.push_back((*itLeft));
                        itLeft++;
                    }

                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (!itRightS->terminal && !itRightS->context) {
                                itRightS++;
                                std::vector<Symbol> rightContextLine;
                                if ( contextSize.second > 0) {
                                    rightContextLine.push_back((*itRightS));
                                    rightContextLine.front().context = true;
                                    if (!rightContext.empty()) {
                                        rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                                        rightContextLine.pop_back();
                                    }
                                }
                                int leftContextId = convertContextToID(0,leftContext);
                                int rightContextID = convertContextToID(1,rightContext);
                                itRightS--;
                                double bInside = p[sizeLeftContext*(*itRightS).id + leftContextId][j][j+k][convertContextToID(1,rightContextLine)];
                                itRightS++;
                                double cInside = p[sizeLeftContext*(*itRightS).id + leftContextId][j+k+1][j+i][rightContextID];
                                p[sizeLeftContext*nonterminal.id+ leftContextId][j][i+j][rightContextID] =
                                        p[sizeLeftContext*nonterminal.id+ leftContextId][j][i+j][rightContextID] + bInside*cInside*(*itRight).second.first;

                                /*std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,rightContext)<<"] = "
                                         << "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,rightContext)<<
                                         "] + IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<j+k<<"]["<< convertContextToID(1,rightContext)<<
                                         "] * IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,leftContext)<<"]["<<j+k+1<<"]["<<j+i<<"]["<< convertContextToID(1,rightContext) << "] * " <<  "ProbR = "
                                         <<p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][j][i+j][convertContextToID(1,rightContext)]<< std::endl;*/
                                //std::cout << p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][j][i+j][convertContextToID(1,rightContext)]  << " + " << bInside << " * " << cInside << " * " << (*itRight).second.first << std::endl;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return p;
}

void Grammar::sampleParseTreeKLVec(
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> &vr,
        Rule r, std::vector<Symbol> w, double ****insideTable, unsigned int i, unsigned int k) {
    unsigned int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    Symbol nonterminal =  Symbol("", 0,false);
    std::vector<Symbol>::iterator itLeft;
    std::vector<Symbol> leftContext;
    for (itLeft = r.left.begin(); itLeft != r.left.end(); itLeft++) {
        if (!itLeft->terminal && !itLeft->context) {
            nonterminal = (*itLeft).clone();
            itLeft++;
            break;
        }
        else {
            leftContext.push_back((*itLeft));
        }
    }
    std::vector<Symbol> rightContext;
    while (itLeft != r.left.end() ) {
        rightContext.push_back((*itLeft));
        itLeft++;
    }
    double sumJBC = 0.0;
    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    std::vector<Symbol> rightContextLine;
                    itSymbol++;
                    if ( contextSize.second > 0) {
                        rightContextLine.push_back((*itSymbol));
                        rightContextLine.front().context = true;
                        if (!rightContext.empty()) {
                            rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                            rightContextLine.pop_back();
                        }
                    }
                    itSymbol--;
                    for (int l = 0; l < jRange; l++) {
                        int unsigned indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                        [i][i + l][convertContextToID(1,rightContextLine)];
                        itSymbol++;
                        double pCjk = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                        [i + l + 1][k][convertContextToID(1,rightContext)];
                        int unsigned indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) /insideTable
                        [contextAmount.first*nonterminal.id+ convertContextToID(0,leftContext)]
                        [i][k][convertContextToID(1,rightContext)];
                        itSymbol--;
                        /*std::cout << "SumJBC = " <<sumJBC << " ";

                        std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first <<
                        " * " << pBij << " ("<< "IT["<<
                                      leftContext.size()*(*itSymbol).id+ convertContextToID(0,leftContext)<<"]["<<i<<"]["<<i+l<<"]["<<
                                      convertContextToID(1,rightContext)<<"])"
                        << " * " << pCjk
                                << " ("<< "IT["<<
                                leftContext.size()*(*itSymbol).id+ convertContextToID(0,leftContext)<<"]["<<i+l+1<<"]["<<k<<"]["<<
                                convertContextToID(1,rightContext)<<"])"
                        << " / " << insideTable
                        [leftContext.size()*nonterminal.id+ convertContextToID(0,leftContext)]
                        [i][k][convertContextToID(1,rightContext)] << " IT["<<
                        leftContext.size()*nonterminal.id+ convertContextToID(0,leftContext)<<"]["<<i<<"]["<<k<<"]["<<
                        convertContextToID(1,rightContext)<<"]"<<std::endl;*/


                    }
                    break;
                }
            }
        }


        int l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //std::cout << "Sensitive: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        unsigned int j = 0;

        delete[]jbc;
        Symbol rightB = Symbol("", 0,false);
        Symbol rightC = Symbol("", 0,false);

        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    if ((*itSymbol).id == l / (nNonTerminals*jRange)) {
                        itSymbol++;
                        if(itSymbol == (*itRight).first.end() || (*itSymbol).context)
                            break;
                        if((*itSymbol).id == (l%(nNonTerminals*jRange)) / jRange) {
                            itSymbol--;
                            rightProduction.first.push_back((*itSymbol));
                            rightB = (*itSymbol).clone();
                            j = l - nNonTerminals*jRange* (*itSymbol).id;
                            itSymbol++;
                            if (itSymbol->context)
                                std::cout << "Error\n";
                            rightProduction.first.push_back((*itSymbol));
                            rightC = (*itSymbol).clone();
                            if ((*itSymbol).id > 0)
                                j = j % (*itSymbol).id;
                            rightProduction.second.first = (*itRight).second.first;
                            rightProduction.second.second = (*itRight).second.second;
                            break;
                        } else {
                            itSymbol--;
                            break;
                        }
                    }
                }
            }
        }

        rules[r.index].ruleFrequence[rightB.id*nNonTerminals + rightC.id].second++;
        rightProduction.first.insert(rightProduction.first.end(), rightContext.begin(), rightContext.end());
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        applyProduction(production);
        std::vector<Symbol> leftContext2;
        std::vector<Symbol> rightContext2;
        getActualContext(leftContext2,rightContext2 );

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();

        unsigned int lc , rc;
        std::uniform_int_distribution<int> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<int> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = ((int) rand()%(leftContext2.size()+1));
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::vector<Symbol> lhs;
        lhs.insert(lhs.end(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightB);
        lhs.insert(lhs.end(), rightContext2.begin(), rightContext2.end());
        Rule nextRule = findRuleByLHS(lhs);
        sampleParseTreeKLVec(vr, nextRule, w, insideTable, i, i+j);

        getActualContext(leftContext2,rightContext2 );

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();


        lc = distIntL(mt);
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();

        lhs.clear();
        lhs.insert(lhs.begin(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightC);
        lhs.insert(lhs.end(), rightContext2.begin(), rightContext2.end());
        sampleParseTreeKLVec(vr, findRuleByLHS(lhs), w, insideTable, i+j+1, k);

    } else {
        std::vector<Symbol> leftContext2;
        std::vector<Symbol> rightContext2;
        getActualContext(leftContext2, rightContext2);

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();
        unsigned int lc , rc;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<int> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<int> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w[i].name;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            std::vector<Symbol>::iterator itSymbol;
            bool flagContext = false;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).name == s && !(*itSymbol).context ) {
                    rightProduction.first.push_back((*itSymbol));
                    rightProduction.second.first = (*itRight).second.first;
                    rightProduction.second.second = (*itRight).second.second;
                    flagContext = true;
                    rules[r.index].ruleFrequence[nNonTerminals*nNonTerminals + (*itSymbol).id].second++;
                } else if (flagContext)
                    rightProduction.first.push_back((*itSymbol));
            }

        }
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        applyProduction(production);
        return;
    }
}

std::pair<double,double> Grammar::perplexity(const std::vector<std::vector<Symbol>>& testData, bool normalized) {
    std::pair<double,double> pXTheta = std::make_pair(1.0, 0.0);
    for (const auto& w: testData) {
        double ***iTable  = CYKProbVec(w);
        double pX = iTable[0][0][w.size()-1] ;
        double pXN = pX;
        double ***pN = CYKProbNVec(w);
        pXN /= pN[0][0][w.size()-1];
        freeInsideTable(pN, w.size());
        freeInsideTable(iTable, w.size());
        pXTheta.first += log(pX);
        pXTheta.second += log(pXN);

    }
    return std::make_pair(exp(-(pXTheta.first/testData.size())), exp(-(pXTheta.second/testData.size())));
}

std::pair<double,double> Grammar::perplexityKL(const std::vector<std::vector<Symbol>>& testData, bool normalized) {
    std::pair<double,double> pXTheta = std::make_pair(1.0, 0.0);
    for (const auto& w: testData) {
        double ****iTable  = CYKProbKLVec(w);
        double pX = iTable[0][0][w.size()-1][0];
        double pXN = pX;
        double ****pN = CYKProbKLNVec(w);
        pXN /= pN[0][0][w.size()-1][0];
        freeInsideTableKL(pN, w.size());
        freeInsideTableKL(iTable, w.size());
        pXTheta.first += log(pX);
        pXTheta.second += log(pXN);

    }
    return std::make_pair(exp(-(pXTheta.first/testData.size())), exp(-(pXTheta.second/testData.size())));
}

double ***Grammar::CYKProbNVec(const std::vector<Symbol>& w) {
    auto ***p = new double**[nonterminals.size()];
    for (int i = 0; i < nonterminals.size(); i++) {
        p[i] = new double *[w.size()];
        for (int j = 0; j < w.size(); j++)
            p[i][j] = new double[w.size()];
    }

    for (int i = 0; i < nonterminals.size(); i++)
        for (int j = 0; j < w.size(); j++)
            for (int k = 0; k < w.size(); k++)
                p[i][j][k] = 0.0;

    std::vector<Rule>::iterator itRule;
    for (int i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol>::iterator itLeft;
            Symbol nonterminal =  Symbol("", 0,false);
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context)
                    nonterminal = (*itLeft).clone();
                break;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                        p[nonterminal.id][i][i] += itRight->second.first;
                        //std::cout<< "IT["<<nonterminal.id<<"]["<<i<<"]["<<i<<"] = " << itRight->second.first << " probR" << std::endl;
                    }
                }
            }

        }
    }
    for (int i = 1; i < w.size(); i++) {
        for (int j = 0; j < w.size()-i; j++) {
            for (int k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol>::iterator itLeft;
                    Symbol nonterminal =  Symbol("", 0,false);
                    for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                        if (!itLeft->terminal && !itLeft->context)
                            nonterminal = (*itLeft).clone();
                        break;
                    }
                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (!itRightS->terminal && !itRightS->context) {
                                double bInside = p[(*itRightS).id][j][j+k];
                                itRightS++;
                                double cInside = p[(*itRightS).id][j+k+1][j+i];
                                p[nonterminal.id][j][i+j] = p[nonterminal.id][j][i+j] + bInside*cInside*(*itRight).second.first;
                                /*std::cout<< "IT["<<nonterminal.id<<"]["<<j<<"]["<<i+j<<"] = "
                                << " IT["<<nonterminal.id<<"]["<<j<<"]["<<i+j<<"] + IT[" <<(*itRightS).id<<"]["<<j<<"]["<<j+k<<"]" <<
                                " * IT[" <<(*itRightS).id<<"]["<<j+k+1<<"]["<<j+i<<"]" << " * " <<  "ProbR = "
                                << p[nonterminal.id][j][i+j] << std::endl;*/
                                //std::cout << p[nonterminal.id][j][i+j]  << " + " << p[(*itRightS).id][j][j+k] << " * " << p[(*itRightS).id][j+k+1][j+i] << " * " << (*itRight).second.first << std::endl;
                                break;
                            }
                        }
                    }

                }
            }
        }
    }

    //printInsideTable(p,w.size());
    return p;
}

double ****Grammar::CYKProbKLNVec(const std::vector<Symbol>& w) {
    int sizeLeftContext = 0;
    int sizeRightContext = 0;
    for (int i = 0; i <= contextSize.first; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();

    }
    contextAmount.first = sizeLeftContext;
    for (int i = 0; i <= contextSize.second; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, nonterminals, i, word, true);
//        rightContexts.insert(rightContexts.end(), permutationsTerminals.begin(), permutationsTerminals.end());
        sizeRightContext += permutationsTerminals.size();
    }
    contextAmount.second = sizeRightContext;
    auto ****p = new double***[nonterminals.size()*sizeLeftContext];
    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++) {
        p[i] = new double **[w.size()];
        for (int j = 0; j < w.size(); j++) {
            p[i][j] = new double *[w.size()];
            for (int k = 0; k < w.size(); k++)
                p[i][j][k] = new double[sizeRightContext];
        }
    }

    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++)
        for (int j = 0; j < w.size(); j++)
            for (int k = 0; k < w.size(); k++)
                for (int l = 0; l < sizeRightContext; l++)
                    p[i][j][k][l] = 0.0;

    std::vector<Rule>::iterator itRule;
    for (int i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol>::iterator itLeft;
            Symbol nonterminal =  Symbol("", 0,false);
            std::vector<Symbol> leftContext;
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context) {
                    nonterminal = (*itLeft).clone();
                    itLeft++;
                    break;
                }
                else {
                    leftContext.push_back((*itLeft));
                }
            }
            std::vector<Symbol> rightContext;
            while (itLeft != (*itRule).left.end() ) {
                rightContext.push_back((*itLeft));
                itLeft++;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                            p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][i][i][convertContextToID(1,rightContext)] += itRight->second.first;
                            //std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<i<<"]["<<i<<"]["<< convertContextToID(1,rightContext)<<"] = ProbR"
                            /* <<(*itRight).second.first << std::endl;*/
                    }
                }
            }

        }
    }
    for (int i = 1; i < w.size(); i++) {
        for (int j = 0; j < w.size()-i; j++) {
            for (int k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol>::iterator itLeft;
                    Symbol nonterminal =  Symbol("", 0,false);
                    std::vector<Symbol> leftContext;
                    for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                        if (!itLeft->terminal && !itLeft->context) {
                            nonterminal = (*itLeft).clone();
                            itLeft++;
                            break;
                        }
                        else {
                            leftContext.push_back((*itLeft));
                        }
                    }
                    std::vector<Symbol> rightContext;
                    while (itLeft != (*itRule).left.end() ) {
                        rightContext.push_back((*itLeft));
                        itLeft++;
                    }

                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (!itRightS->terminal && !itRightS->context) {
                                itRightS++;
                                std::vector<Symbol> rightContextLine;
                                if ( contextSize.second > 0) {
                                    rightContextLine.push_back((*itRightS));
                                    rightContextLine.front().context = true;
                                    if (!rightContext.empty()) {
                                        rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                                        rightContext.pop_back();
                                    }
                                }
                                itRightS--;
                                int leftContextID = convertContextToID(0,leftContext);
                                int rightContextID = convertContextToID(1,rightContext);
                                double bInside = p[sizeLeftContext*(*itRightS).id + leftContextID][j][j+k][convertContextToID(1,rightContextLine)];
                                itRightS++;
                                double cInside = p[sizeLeftContext*(*itRightS).id + leftContextID][j+k+1][j+i][rightContextID];
                                p[sizeLeftContext*nonterminal.id+ leftContextID][j][i+j][rightContextID] =
                                        p[sizeLeftContext*nonterminal.id+ leftContextID][j][i+j][rightContextID] + bInside*cInside*(*itRight).second.first;

                                /*std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,rightContext)<<"] = "
                                         << "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,rightContext)<<
                                         "] + IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,leftContext)<<"]["<<j<<"]["<<j+k<<"]["<< convertContextToID(1,rightContext)<<
                                         "] * IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,leftContext)<<"]["<<j+k+1<<"]["<<j+i<<"]["<< convertContextToID(1,rightContext) << "] * " <<  "ProbR = "
                                         <<p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][j][i+j][convertContextToID(1,rightContext)]<< std::endl;*/
                                //std::cout << p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][j][i+j][convertContextToID(1,rightContext)]  << " + " << bInside << " * " << cInside << " * " << (*itRight).second.first << std::endl;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return p;
}

Grammar::~Grammar() = default;

void Grammar::gibbsSamplingPCFG(int iterations) {
    //srand((unsigned) time(nullptr));
    int sentencesLength = words.size();
    for (int j = 0; j < iterations; j++) {
        parseTreesVec.clear();
        for (int i = 0; i< sentencesLength; i++) {
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
            double ***insideTable;
            insideTable = CYKProbVec(words[i]);
            sampleParseTreeVec(vProds, rules[0], words[i], insideTable, 0, words[i].size()-1);
            freeInsideTable(insideTable, words[i].size());
            std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTree;
            pairTree.first = words[i];
            pairTree.second = vProds;
            parseTreesVec.push_back(pairTree);
        }
        std::vector<Symbol> a;
        calculateNewThetaVecOpt(a, -1);
        updateParseTressTheta();
        if (j%(iterations/10) == 0) {
            std::pair<double,  double> pMpP = perplexity(words, true);
            std::cout << "iteration: " << j <<" - PerplexityM: "<< pMpP.first << " - PerplexityP: "<< pMpP.second << std::endl;
        }
    }
}

void Grammar::gibbsSamplingPCSG(int iterations) {
    auto totalTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto insideTTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto iterationTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto newThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto updateThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto perplexityTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    //srand((unsigned) time(nullptr));
    int sentencesLength = words.size();
    auto startMetropolis = std::chrono::system_clock::now();
    for (int j = 0; j < iterations; j++) {
        parseTreesVec.clear();
        for (int i = 0; i< sentencesLength; i++) {
            auto startIt = std::chrono::system_clock::now();
            actualProduction.clear();
            actualProduction.push_back(nonterminals[0]);
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
            auto startInsideTable = std::chrono::system_clock::now();
            double ****insideTableKL;
            insideTableKL = CYKProbKLVec(words[i]);
            insideTTime += std::chrono::system_clock::now() - startInsideTable;
            sampleParseTreeKLVec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size()-1);
            freeInsideTableKL(insideTableKL, words[i].size());
            std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTreeVec;
            pairTreeVec.first = words[i];
            pairTreeVec.second = vProds;
            parseTreesVec.push_back(pairTreeVec);

            iterationTime = std::chrono::system_clock::now() - startIt;
            //std::cout << "   insideT: " << insideTTime.count() << " iterationT: " << iterationTime.count() << std::endl;
        }
        std::vector<Symbol> a;
        auto startNewTheta = std::chrono::system_clock::now();
        calculateNewThetaVecOpt(a, -1);
        newThetaTime += std::chrono::system_clock::now() - startNewTheta;
        auto startUpdateTheta = std::chrono::system_clock::now();
        updateParseTressTheta();
        updateThetaTime += std::chrono::system_clock::now() - startUpdateTheta;

        //std::cout << "   newTheta: " << newThetaTime.count() << " updateTheta: " << updateThetaTime.count() << std::endl;

        if (j%(iterations/5) == 0) {
            //auto startPerplexityTime = std::chrono::system_clock::now();
            //std::pair<double,  double> pMpP = perplexityKL(words, true);
            //perplexityTime += std::chrono::system_clock::now() - startPerplexityTime;
            //std::cout << "      iteration: " << j <<" - Perplexity: "<< pMpP.first << " - PerplexityN: "<< pMpP.second <<std::endl;
            std::cout << "      iteration: " << j  << std::endl;
        }
    }
}

void Grammar::freeInsideTable(double ***p, int wSize) const {
    for (int i = 0; i < nonterminals.size(); i++) {
        for (int j = 0; j < wSize; j++)
            delete[] p[i][j];
        delete[] p[i];
    }
    delete[] p;
}

void Grammar::freeInsideTableKL(double ****p, int wSize) const {
    for (int i = 0; i < nonterminals.size()*contextAmount.first; i++) {
        for (int j = 0; j < wSize; j++) {
            for (int k = 0; k < wSize; k++)
                delete[] p[i][j][k];
            delete[] p[i][j];
        }
        delete[] p[i];
    }
    delete p;
}

void Grammar::sampleParseTreeKLVecOpt(
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> &vr,
        Rule r, std::vector<Symbol> w, double ****insideTable, unsigned int i, unsigned int k) {
    unsigned int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    Symbol nonterminal =  Symbol("", 0,false);
    std::vector<Symbol>::iterator itLeft;
    std::vector<Symbol> leftContext;
    for (itLeft = r.left.begin(); itLeft != r.left.end(); itLeft++) {
        if (!itLeft->terminal && !itLeft->context) {
            nonterminal = (*itLeft).clone();
            itLeft++;
            break;
        }
        else {
            leftContext.push_back((*itLeft));
        }
    }
    std::vector<Symbol> rightContext;
    while (itLeft != r.left.end() ) {
        rightContext.push_back((*itLeft));
        itLeft++;
    }
    double sumJBC = 0.0;
    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    std::vector<Symbol> rightContextLine;
                    itSymbol++;
                    if ( contextSize.second > 0) {
                        rightContextLine.push_back((*itSymbol));
                        rightContextLine.front().context = true;
                        if (!rightContext.empty()) {
                            rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                            rightContext.pop_back();
                        }
                    }
                    itSymbol--;
                    for (int l = 0; l < jRange; l++) {
                        unsigned int indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                        [i][i + l][convertContextToID(1,rightContextLine)];
                        itSymbol++;
                        double pCjk = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                        [i + l + 1][k][convertContextToID(1,rightContext)];
                        unsigned int indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) /insideTable
                        [contextAmount.first*nonterminal.id+ convertContextToID(0,leftContext)]
                        [i][k][convertContextToID(1,rightContext)];
                        itSymbol--;
                        /*std::cout << "SumJBC = " <<sumJBC << " ";

                        std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first <<
                        " * " << pBij << " ("<< "IT["<<
                                      leftContext.size()*(*itSymbol).id+ convertContextToID(0,leftContext)<<"]["<<i<<"]["<<i+l<<"]["<<
                                      convertContextToID(1,rightContext)<<"])"
                        << " * " << pCjk
                                << " ("<< "IT["<<
                                leftContext.size()*(*itSymbol).id+ convertContextToID(0,leftContext)<<"]["<<i+l+1<<"]["<<k<<"]["<<
                                convertContextToID(1,rightContext)<<"])"
                        << " / " << insideTable
                        [leftContext.size()*nonterminal.id+ convertContextToID(0,leftContext)]
                        [i][k][convertContextToID(1,rightContext)] << " IT["<<
                        leftContext.size()*nonterminal.id+ convertContextToID(0,leftContext)<<"]["<<i<<"]["<<k<<"]["<<
                        convertContextToID(1,rightContext)<<"]"<<std::endl;*/


                    }
                    break;
                }
            }
        }


        int l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //std::cout << "Sensitive: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        unsigned int j = 0;

        delete[]jbc;
        Symbol rightB = Symbol("", 0,false);
        Symbol rightC = Symbol("", 0,false);

        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        rightProduction = r.getRightSidebyId(l / (nNonTerminals*jRange), (l%(nNonTerminals*jRange)) / jRange, false, nNonTerminals);
        std::vector<Symbol>::iterator itSymbol;
        for (itSymbol = rightProduction.first.begin(); itSymbol != rightProduction.first.end(); itSymbol++) {
            if (!(*itSymbol).terminal && !(*itSymbol).context) {
                rightB = (*itSymbol).clone();
                j = l - nNonTerminals*jRange* (*itSymbol).id;
                itSymbol++;
                rightC = (*itSymbol).clone();
                if ((*itSymbol).id > 0)
                    j = j % (*itSymbol).id;
                break;
            }
        }
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);

        vr.push_back(production);
        applyProduction(production);
        std::vector<Symbol> leftContext2;
        std::vector<Symbol> rightContext2;
        getActualContext(leftContext2,rightContext2 );

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();

        unsigned int lc , rc;
        std::uniform_int_distribution<int> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<int> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = ((int) rand()%(leftContext2.size()+1));
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::vector<Symbol> lhs;
        lhs.insert(lhs.end(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightB);
        lhs.insert(lhs.end(), rightContext2.begin(), rightContext2.end());
        Rule nextRule = findRuleByLHS(lhs);
        sampleParseTreeKLVecOpt(vr, nextRule, w, insideTable, i, i+j);

        getActualContext(leftContext2,rightContext2 );

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();


        lc = distIntL(mt);
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();

        lhs.clear();
        lhs.insert(lhs.begin(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightC);
        lhs.insert(lhs.begin(), rightContext2.begin(), rightContext2.end());
        sampleParseTreeKLVecOpt(vr, findRuleByLHS(lhs), w, insideTable, i+j+1, k);

    } else {
        std::vector<Symbol> leftContext2;
        std::vector<Symbol> rightContext2;
        getActualContext(leftContext2, rightContext2);

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();
        unsigned int lc , rc;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<int> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<int> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w[i].name;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            std::vector<Symbol>::iterator itSymbol;
            bool flagContext = false;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).name == s && !(*itSymbol).context ) {
                    rightProduction.first.push_back((*itSymbol));
                    rightProduction.second.first = (*itRight).second.first;
                    rightProduction.second.second = (*itRight).second.second;
                    flagContext = true;
                } else if (flagContext)
                    rightProduction.first.push_back((*itSymbol));
            }

        }
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        applyProduction(production);
        return;
    }
}

double ****Grammar::CYKProbKLVecOpt(std::vector<Symbol> w) {
    int sizeLeftContext = 0;
    int sizeRightContext = 0;
    for (int i = 0; i <= contextSize.first; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();

    }
    contextAmount.first = sizeLeftContext;
    for (int i = 0; i <= contextSize.second; i++) {
        std::vector<Symbol> word;
        std::vector<std::vector<Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, nonterminals, i, word, true);
//        rightContexts.insert(rightContexts.end(), permutationsTerminals.begin(), permutationsTerminals.end());
        sizeRightContext += permutationsTerminals.size();
    }
    contextAmount.second = sizeRightContext;
    auto ****p = new double***[nonterminals.size()*sizeLeftContext];
    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++) {
        p[i] = new double **[w.size()];
        for (int j = 0; j < w.size(); j++) {
            p[i][j] = new double *[w.size()];
            for (int k = 0; k < w.size(); k++)
                p[i][j][k] = new double[sizeRightContext];
        }
    }

    for (int i = 0; i < nonterminals.size()*sizeLeftContext; i++)
        for (int j = 0; j < w.size(); j++)
            for (int k = 0; k < w.size(); k++)
                for (int l = 0; l < sizeRightContext; l++)
                    p[i][j][k][l] = 0.0;
    std::vector<Rule>::iterator itRule;
    for (int i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol>::iterator itLeft;
            Symbol nonterminal =  Symbol("", 0,false);
            std::vector<Symbol> leftContext;
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context) {
                    nonterminal = (*itLeft).clone();
                    itLeft++;
                    break;
                }
                else {
                    leftContext.push_back((*itLeft));
                }
            }
            std::vector<Symbol> rightContext;
            while (itLeft != (*itRule).left.end() ) {
                rightContext.push_back((*itLeft));
                itLeft++;
            }

            std::pair<std::vector<Symbol>,std::pair<double, double>> rightHandSide;
            rightHandSide = (*itRule).getRightSidebyId(w[i].id, 0, true, nNonTerminals);
            std::vector<Symbol>::iterator itRightS;
            for (itRightS = rightHandSide.first.begin(); itRightS != rightHandSide.first.end(); itRightS++) {
                if (itRightS->terminal && !itRightS->context) {
                    if(itRightS->name == w[i].name) {
                        p[sizeLeftContext*nonterminal.id+ convertContextToID(0,leftContext)][i][i][convertContextToID(1,rightContext)] = rightHandSide.second.first;
                        /*std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,leftContext)<<"]["<<i<<"]["<<i<<"]["<< convertContextToID(1,rightContext)<<"] = ProbR"
                        << " " << (*itRight).second.first << std::endl;*/
                        break;
                    }
                }
            }

        }
    }
    for (int i = 1; i < w.size(); i++) {
        for (int j = 0; j < w.size()-i; j++) {
            for (int k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                            if (!(*itRight).first[(*itRule).index1stNonContext].terminal &&
                                !(*itRight).first[(*itRule).index1stNonContext].context) {
                                std::vector<Symbol> rightContextLine;
                                if ( contextSize.second > 0) {
                                    rightContextLine.push_back((*itRight).first[(*itRule).index1stNonContext+1]);
                                    rightContextLine.front().context = true;
                                    if (!(*itRule).rightContext.empty()) {
                                        rightContextLine.insert(rightContextLine.end(), (*itRule).rightContext.begin(), (*itRule).rightContext.end());
                                        rightContextLine.pop_back();
                                    }
                                }
                                int leftContextID = convertContextToID(0,(*itRule).leftContext);
                                int rightContextID = convertContextToID(1,(*itRule).rightContext);
                                double bInside = p[sizeLeftContext*(*itRight).first[(*itRule).index1stNonContext].id + leftContextID][j][j+k][convertContextToID(1,rightContextLine)];
                                double cInside = p[sizeLeftContext*(*itRight).first[(*itRule).index1stNonContext+1].id + leftContextID][j+k+1][j+i][rightContextID];

                                p[sizeLeftContext*(*itRight).first[(*itRule).index1stNonContext].id+ leftContextID][j][i+j][rightContextID] =
                                        p[sizeLeftContext*(*itRight).first[(*itRule).index1stNonContext].id+ leftContextID][j][i+j][rightContextID] + bInside*cInside*(*itRight).second.first;



/*                                for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                                    std::vector<Symbol>::iterator itRightS;
                                    std::vector<Symbol> rightContextLine;
                                    if ( contextSize.second > 0) {
                                        rightContextLine.push_back((*itRight).first[(*itRule).index1stNonContext+1]);
                                        rightContextLine.front().context = true;
                                        if (!(*itRule).rightContext.empty()) {
                                            rightContextLine.insert(rightContextLine.end(), (*itRule).rightContext.begin(), (*itRule).rightContext.end());
                                            rightContextLine.pop_back();
                                        }
                                    }
                                    int leftContextID = convertContextToID(0,(*itRule).leftContext);
                                    int rightContextID = convertContextToID(1,(*itRule).rightContext);
                                    double bInside = p[sizeLeftContext*((*itRight).first[(*itRule).index1stNonContext]).id + leftContextID][j][j+k][convertContextToID(1,rightContextLine)];
                                    double cInside = p[sizeLeftContext*((*itRight).first[(*itRule).index1stNonContext+1]).id + leftContextID][j+k+1][j+i][rightContextID];
                                    //nonterminal  é igual a (*itRight).first[(*itRule).index1stNonContext]
                                    p[sizeLeftContext*(*itRight).first[(*itRule).index1stNonContext].id+ leftContextID][j][i+j][rightContextID] =
                                            p[sizeLeftContext*(*itRight).first[(*itRule).index1stNonContext].id+ leftContextID][j][i+j][rightContextID] + bInside*cInside*(*itRight).second.first;
                                    break;
                                }*/

                        }
                    }
                }
            }
        }
    }
    //printInsideTableKL(p,w.size());
    return p;
}


void Grammar::generateRulesRegular() {
    for (const auto& n: nonterminals) {
        std::vector<Symbol> left;
        left.push_back(n);
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>> right;
        for (const auto& t: terminals) {
            std::pair<std::vector<Symbol>,std::pair<double, double>> rightHS;
            rightHS.first.push_back(t);
            for (const auto& nr: nonterminals) {
                rightHS.first.push_back(nr);
                rightHS.second.first = 1.0/(nNonTerminals*nTerminals+1);
                rightHS.second.second = ALFA;
                right.push_back(rightHS);
                rightHS.first.pop_back();
            }
        }
        std::pair<std::vector<Symbol>, std::pair<double, double>> rightHS;
        rightHS.second.first = 1.0/(nNonTerminals*nTerminals+1);
        rightHS.second.second = ALFA;
        right.push_back(rightHS);
        Rule r = Rule(left, right);
        rules.push_back(r);
    }
}

void Grammar::baumWelch(int iteration) {
    std::vector<double> countNd;
    std::vector<double> countNl;
    vector<vector<vector<double>>> countNT;
    for (int i = 0; i < iteration; i++) {
        cout << "Iteration: "<<i << endl;
        baumWelchExpectation(countNd, countNl, countNT);
        baumWelchMaximization(countNd, countNl, countNT);
    }

}


void Grammar::baumWelchExpectation(std::vector<double> &countNd, std::vector<double> &countNl, std::vector<std::vector<std::vector<double>>> &countNT) {
    countNd.clear();
    countNl.clear();
    countNT.clear();
    for (const auto& n: nonterminals) {
        countNd.push_back(0);
        countNl.push_back(0);
        countNT.emplace_back();
        for (const auto& t: terminals) {
            countNT[n.id].push_back(vector<double>());
            for (auto nl: nonterminals) {
                countNT[n.id][t.id].push_back(0.0);
            }
        }
    }
    for (auto x: words) {
        vector<vector<double>> F, B;
        for (int i = 0; i <= x.size(); i++) {
            F.emplace_back();
            B.emplace_back();
            for (auto n: nonterminals) {
                F[i].push_back(0.0);
                B[i].push_back(0.0);
            }
        }
        //To DO: contar quantas vezes x acontece ou noã precisa?
        for (const auto& n: nonterminals) {
            if (n.id == 0) F[0][n.id] = 1; else F[0][n.id] = 0;
            B[x.size()][n.id] = rules[n.id].right[nNonTerminals * nTerminals].second.first;
            //cout << "F[" << "0" << "][" << n.id << "] = " << F[0][n.id] << "\t\t\t\t\t\t\t\t";
            //cout << "B[" << x.size() << "][" << n.id << "] = "
            //     << rules[n.id].right[nNonTerminals * nTerminals].second.first << endl;
        }
        for (int i = 1; i<= x.size(); i++) {
            for (const auto& n: nonterminals) {
                F[i][n.id] = 0;
                B[x.size()-i][n.id]=0;
                for (const auto& nl: nonterminals) {
                    F[i][n.id] += F[i-1][nl.id]*rules[nl.id].right[x[i-1].id*nTerminals+n.id].second.first;
                //    cout << "F["<<i<<"]"<<"["<<n.id<<"]"<< " += F["<<i-1<<"]"<<"["<<nl.id<<"] * " << "P["<<nl.id<<"]["<<x[i-1].id<<"]["<<n.id<<"] = "
                //    << F[i-1][nl.id] << "*" << rules[nl.id].right[x[i-1].id*nTerminals+n.id].second.first << " = " << F[i][n.id]<< "\t\t\t\t";
                    B[x.size()-i][n.id] +=  B[x.size()-i+1][nl.id]* rules[n.id].right[x[x.size()-i].id*nTerminals+nl.id].second.first;
                //    cout << "B["<<x.size()-i<<"]"<<"["<<n.id<<"]"<< " += " << "B["<<x.size()-i+1<<"]"<<"["<<nl.id<<"]" << " * P["<<n.id<<"]["<<x[i-1].id<<"]["<<nl.id<<"] = "
                //    << B[x.size()-i+1][nl.id]<< "*" << rules[n.id].right[x[x.size()-i].id*nTerminals+nl.id].second.first << " = " << B[x.size()-i][n.id] <<endl;
                }
            }
        }
        double total = 0.0;
        double totalb = B[0][0];
        for (const auto& n: nonterminals) {
            total += F[x.size()][n.id]*rules[n.id].right[nNonTerminals*nTerminals].second.first;
            //cout << " total = " << "total + " <<  F[x.size()][n.id] << "*" << rules[n.id].right[nNonTerminals*nTerminals].second.first << endl;
        }
        //cout << "Total = " << total << endl;
        //cout << "TotalB = " << totalb <<endl;
        for (int i = 1; i<= x.size(); i++) {
            for (const auto& n: nonterminals) {
                for (const auto& nl: nonterminals) {
                    double val =(F[i-1][n.id]*  rules[n.id].right[x[i-1].id*nTerminals+nl.id].second.first* B[i][nl.id])/total;
                    //cout << "val = F[" <<i-1<<"]["<<n.id<<"]* P["<<n.id<<"]["<<x[i-1].id<<"]["<<nl.id<<"]*"<<"B["<<i<<"]["<<nl.id<<"] = " <<
                    //F[0][0] << "*" << rules[n.id].right[x[i-1].id*nTerminals+nl.id].second.first <<"*" << B[i][nl.id]<< " = " << val << " / ";
                    countNT[n.id][x[i-1].id][nl.id] += val;
                    //cout << "countNT["<<n.id<<"]["<<x[i-1].id<<"]["<<nl.id<<"] += val  = " <<countNT[n.id][x[i-1].id][nl.id] << " / ";
                    countNl[nl.id] += val;
                    //cout << "countNl["<<nl.id<<"] += "<< countNl[nl.id] << endl;
                }

            }
        }
        for (const auto& n: nonterminals) {
            countNd[n.id] += (F[x.size()][n.id]* rules[n.id].right[nNonTerminals*nTerminals].second.first)/total;
            //cout << "countNd["<<n.id<<"] += F["<<x.size()<<"]["<<n.id<<"]"<<"*"<< "P["<<n.id<<"][F]["<<nNonTerminals*nTerminals<<"] = "
            //<< F[x.size()][n.id] << "*" <<rules[n.id].right[nNonTerminals*nTerminals].second.first << " / " << total << " = " << countNd[n.id] << endl;
            if (n.id == 0) countNl[n.id] += B[0][n.id]/total;
        }
    }
}

void Grammar::baumWelchMaximization(vector<double> &countNdR, vector<double> &countNlR,
                               vector<std::vector<std::vector<double>>> &countNTR) {
    vector<double> val;
    for (auto n: nonterminals) {
        val.push_back(0.0);
    }
    for (const auto& n: nonterminals) {
        rules[n.id].right[nNonTerminals*nTerminals].second.first = (1.0*countNdR[n.id])/(1.0*countNlR[n.id]);
        //cout << "rulesF["<<n.id<<"]["<<nNonTerminals*nTerminals<<"] = countNdR["<<n.id<<"]/countNlR["<<n.id<<"] = " << countNdR[n.id] <<"/"<< countNlR[n.id] << endl;
        double total = 0.0;
        for (const auto& t: terminals) {
            for (const auto& nl: nonterminals) {
                total += 1.0 * countNTR[n.id][t.id][nl.id];
                rules[n.id].right[t.id*nTerminals+nl.id].second.first = (1.0*countNTR[n.id][t.id][nl.id])/(1.0*countNlR[n.id]);
                //cout << "rulesF["<<n.id<<"]["<<t.id*nTerminals+nl.id<<"] = countNTR["<<n.id<<"]["<<t.id<<"]["<<nl.id<<"]/countNlR["<<n.id<<"] = "
                //<< countNTR[n.id][t.id][nl.id] <<"/"<< countNlR[n.id] << endl;
            }
        }
        //cout << "cNdR: " << countNdR[n.id] << " cNTR: " << total
        //<< " sum: " << countNdR[n.id]+ total<< " cNlR: "<< countNlR[n.id] << endl;
    }
}

void Grammar::insideOutside(int iterations) {
    for (int i = 0; i < iterations; i++) {
        vector<vector<vector<double>>>newProbsNTn;
        vector<vector<vector<double>>>newProbsNTd;
        vector<vector<double>>newProbsTn;
        for (int i2 = 0 ; i2 <nNonTerminals; i2++) {
            newProbsNTn.emplace_back();
            newProbsNTd.emplace_back();
            for (int j = 0; j < nNonTerminals; j++) {
                newProbsNTn[i2].push_back(vector<double>());
                newProbsNTd[i2].push_back(vector<double>());
                for (int k = 0; k < nNonTerminals; k++) {
                    newProbsNTn[i2][j].push_back(0.0);
                    newProbsNTd[i2][j].push_back(0.0);
                }
            }
        }
        for (int i2 = 0 ; i2 < nNonTerminals; i2++) {
            newProbsTn.emplace_back();
            for (int m = 0; m < nTerminals; m++) {
                newProbsTn[i2].push_back(0.0);
            }
        }

        cout << "Iteration: "<<i << endl;
        for (auto w: words) {
            double *** insideTable;
            insideTable = CYKProbVec(w);
            //printInsideTable(insideTable, w.size());
            std::vector<std::vector<std::vector<double>>> oT = outsideTable(w, insideTable);
            //printOutsideTable(oT);
            vector<vector<vector<vector<vector<double>>>>> W;
            vector<vector<vector<double>>> V;
            calculateWV(W,V, w, insideTable, oT);

            for (int s = 0 ; s <w.size()-1; s++) {
                for (int t = s+1; t <w.size(); t++) {
                    for (int i2 = 0; i2 < nNonTerminals; i2++) {
                        for (int j = 0; j < nNonTerminals; j++) {
                            for (int k = 0; k < nNonTerminals; k++) {
                                    newProbsNTn[i2][j][k] += W[s][t][i2][j][k];
                                    //cout << "nPNTn["<<i2<<"]["<<j<<"]["<<k<<"] = "<< newProbsNTn[i2][j][k] << endl;
                            }
                        }
                    }
                }
            }

            for (int i2 = 0; i2 < nNonTerminals; i2++) {
                for (int j = 0; j < nNonTerminals; j++) {
                    for (int k = 0; k < nNonTerminals; k++) {
                        for (int s = 0 ; s <w.size(); s++) {
                            for (int t = s; t <w.size(); t++) {
                                    newProbsNTd[i2][j][k] += V[s][t][i2];
                                    //cout << "nPNTd["<<i2<<"]["<<j<<"]["<<k<<"] = "<< newProbsNTd[i2][j][k] << endl;
                            }
                        }
                    }
                }
            }

            for (int t = 0; t <w.size(); t++) {
                for (int i2 = 0; i2 < nNonTerminals; i2++) {
                    for (int m = 0; m < nTerminals; m++) {
                        if(terminals[m].id == w[t].id) {
                            newProbsTn[i2][m] += V[t][t][i2];
                            //cout << "nPT["<<i2<<"]["<<m<<"] = "<< newProbsTn[i2][m] << endl;
                        }
                    }
                }
            }
        }

        vector<Rule>::iterator itRule;
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                if (!(*itRight).first[0].terminal) {
                    (*itRight).second.first = newProbsNTn[(*itRule).left[0].id][(*itRight).first[0].id][(*itRight).first[1].id]/
                                         newProbsNTd[(*itRule).left[0].id][(*itRight).first[0].id][(*itRight).first[1].id];
                }
                else {
                    (*itRight).second.first = newProbsTn[(*itRule).left[0].id][(*itRight).first[0].id]/
                                         newProbsNTd[(*itRule).left[0].id][0][0];
                    //cout << "P"<<(*itRule).left[0].name<<"->"<< (*itRight).first[0].name<< " = " << (*itRight).second.first << endl;
                }
            }
        }
        printGrammar();
    }
}

void Grammar::calculateWV(vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &W,
                          vector<std::vector<std::vector<double>>> &V, const vector<Symbol>& w,double***insideTable, vector<vector<vector<double>>>oT) {

    for (int s = 0 ; s <w.size(); s++) {
        W.emplace_back();
        for (int t = 0; t <w.size(); t++) {
            W[s].push_back(vector<vector<vector<double>>>());
            for (int i = 0; i < nNonTerminals; i++) {
                W[s][t].push_back(vector<vector<double>>());
                for (int j = 0; j < nNonTerminals; j++) {
                    W[s][t][i].push_back(vector<double>());
                    for (int k = 0; k < nNonTerminals; k++) {
                        W[s][t][i][j].push_back(0.0);
                    }
                }
            }
        }
    }

    for (int s = 0 ; s <w.size(); s++) {
        V.emplace_back();
        for (int t = 0; t < w.size(); t++) {
            V[s].push_back(vector<double>());
            for (int i = 0; i < nNonTerminals; i++) {
                V[s][t].push_back(0.0);
            }
        }
    }

    for (int s = 0 ; s <w.size(); s++) {
        for (int t = 0; t <w.size(); t++) {
            for (int i = 0; i < nNonTerminals; i++) {
                for (int j = 0; j < nNonTerminals; j++) {
                    for (int k = 0; k < nNonTerminals; k++) {
                        for (int r = s; r < t; r++) {
                            W[s][t][i][j][k] += (rules[i].right[k+ nNonTerminals * j].second.first *
                                                insideTable[j][s][r] * insideTable[k][r+1][t] * oT[i][s][t]) / insideTable[0][0][w.size()-1];
                            //cout<< "W["<<s<<"]["<<t<<"]["<<i<<"]["<<j<<"]["<<k<<"] = " << W[s][t][i][j][k] << endl;
                        }
                    }
                }
            }
        }
    }

    for (int s = 0 ; s <w.size(); s++) {
        for (int t = 0; t < w.size(); t++) {
            for (int i = 0; i < nNonTerminals; i++) {
                V[s][t][i] = (insideTable[i][s][t]*oT[i][s][t])/insideTable[0][0][w.size()-1];
                //cout<< "V["<<s<<"]["<<t<<"]["<<i<<"] = " << V[s][t][i] << endl;
            }
        }
    }
}

void Grammar::genFPTA() {
    Symbol nI = Symbol("NTI", 0, false, false);
    nonterminals.push_back(nI);
    int nNT = 1;
    vector<Symbol> vr;
    vr.push_back(nI);
    Rule r = Rule(vr, std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>());
    rules.push_back(r);
    for (const auto& w: words) {
        string s;
        for (auto & i : w) {
            s += i.name;
            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            bool existRule = false;
            for (itRight = rules[nI.id].right.begin(); itRight != rules[nI.id].right.end()-1; itRight++) {
                if ((*itRight).first[1].name == "NT"+s) {
                    (*itRight).second.first += 1.0;
                    existRule = true;
                    nI = (*itRight).first[1];
                    break;
                }
            }
            if (!existRule) {
                Symbol nt = Symbol("NT"+s, nNT, false, false);
                nonterminals.push_back(nt);
                nNT++;
                std::pair<std::vector<Symbol>,std::pair<double, double>> right;
                right.first.push_back(i);
                right.first.push_back(nt);
                right.second.first = 1.0;
                rules[nI.id].right.insert(rules[nI.id].right.end()-1, right);
                //rules[nI.id].right.push_back(right);
                vr.clear();
                vr.push_back(nt);
                r = Rule(vr, std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>());
                rules.push_back(r);
                nI = nt;
            }

        }
        bool existRule = false;
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = rules[nI.id].right.begin(); itRight != rules[nI.id].right.end(); itRight++) {
            if ((*itRight).first[0].name.empty()) {
                (*itRight).second.first += 1.0;
                existRule = true;
            }
        }
        if (!existRule) {
            Symbol nt = Symbol("", -1, false, false);
            nonterminals.push_back(nt);
            std::pair<std::vector<Symbol>,std::pair<double, double>> right;
            right.first.push_back(nt);
            right.second.first = 1.0;
            rules[nI.id].right.push_back(right);

        }
        nI = nonterminals[0];
    }
    nNonTerminals = nNT;
}

void Grammar::ALERGIA(double alpha) {
    genFPTA();
    printGrammar();
    vector<Symbol> red;
    red.push_back(nonterminals[0]);
    vector<Symbol> blue;
    for (auto s: rules[0].right) {
        if (s.first[0].name.empty())
            break;
        blue.push_back(s.first[1]);
    }
    int t0 = 30;


    auto itBlue = blue.begin();
    double freqCBlue = rules[(*itBlue).id].freq();
    while (freqCBlue >= t0) {
        bool existRed = false;
        for (const auto& cRed: red) {
            if (compatibleALERGIA(cRed, (*itBlue), alpha)) {
                cout << "merge " << cRed.name << " and " << (*itBlue).name << endl;
                stochasticMerge(cRed, (*itBlue));
                existRed = true;
                break;
            }
        }
        if (!existRed) {
            red.push_back((*itBlue));
            for (auto right: rules[(*itBlue).id].right) {
                if (!(right.first[0].name.empty()))
                    blue.push_back(right.first[1]);
            }
            itBlue = blue.begin();
        }
        blue.erase(itBlue);
        itBlue = blue.begin();
        if (itBlue == blue.end())
            break;
        freqCBlue = rules[(*itBlue).id].freq();
        printGrammar();

    }
    removeUnusedNT();
    normalizeProbs();
}

bool Grammar::compatibleALERGIA(const Symbol& a, const Symbol& b, double alpha) {
    bool correct = true;
    if (!testALERGIA((*(rules[a.id].right.end()-1)).second.first, rules[a.id].freq(),(*(rules[b.id].right.end()-1)).second.first, rules[b.id].freq(), alpha))
        correct = false;
    for (const auto& t: terminals) {
        double dFreqa = 0.0;
        double dFreqb = 0.0;

        for (auto nt: rules[a.id].right) {
            if (nt.first[0].id == t.id) {
                dFreqa =nt.second.first;
            }
        }
        for (auto nt: rules[b.id].right) {
            if (nt.first[0].id == t.id) {
                dFreqa =nt.second.first;
            }
        }
        if (!testALERGIA(dFreqa, rules[a.id].freq(), dFreqb, rules[b.id].freq(), alpha))
            correct = false;
    }
    return correct;
}

bool Grammar::testALERGIA(double probFa, double freqa, double probFb, double freqb, double alpha) {
    double y = abs(probFa/freqa - probFb/freqb);
    return (y < ( (sqrt(1/freqa) + sqrt(1/freqb)) * sqrt(0.5 * log(2/alpha)) ));
}

void Grammar::stochasticMerge(const Symbol& a, const Symbol& b) {
    vector<Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule < rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            if ((*itRight).first[1].id == b.id) {
                stochasticFold(a,b);
                std::pair<std::vector<Symbol>,std::pair<double, double>> newRight;
                int n = (*itRight).second.first;
                newRight.first.push_back((*itRight).first[0]);
                newRight.first.push_back(a);
                newRight.second.first = (*itRight).second.first;
                (*itRight).second.first = 0;
                (*itRule).right.insert((*itRule).right.end()-1, newRight);
                return;
            }
        }
    }
}

void Grammar::stochasticFold(const Symbol& a, const Symbol& b) {
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = rules[b.id].right.begin(); itRight != rules[b.id].right.end()-1; itRight++) {
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRighta;
        for (itRighta = rules[a.id].right.begin(); itRighta != rules[a.id].right.end(); itRighta++) {
            if (itRighta->first[0].id == itRight->first[0].id) {
                if (itRighta->second.first >=1.0)
                    break;
            }
        }
        if (itRighta != rules[a.id].right.end()) {
            if (itRighta->second.first >= 1.0) {
                stochasticFold(itRighta->first[1], itRight->first[1]);
                (*itRighta).second.first += (*itRight).second.first;
            }
        }
        else {
            itRighta->first[1] = Symbol(itRight->first[1]);
            (*itRighta).second.first = (*itRight).second.first;

        }
    }
    (*(rules[a.id].right.end()-1)).second.first += (*(rules[b.id].right.end()-1)).second.first;

}

void Grammar::removeUnusedNT() {
    vector<Symbol> notUsed;
    vector<Symbol> used;
    used.push_back(rules[0].left[0]);
    for (int i = 0; i < used.size(); i++) {
        Rule r = rules[used[i].id];
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end()-1; itRight++) {
            if ((*itRight).second.first <1.0) {
                notUsed.push_back((*itRight).first[1]);
                Symbol aux = (*itRight).first[1];
                recursiveInsertUnused(notUsed, aux);

            }
            else {
                bool isUsed = false;
                for (const auto& u: used) {
                    if(itRight->first[1].id == u.id) {
                        isUsed = true;
                        break;
                    }
                }
                if (!isUsed)
                    used.push_back((*itRight).first[1]);
            }
        }
    }
    for (const auto& nu: notUsed) {
        vector<Rule>::iterator itRule;
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            if ((*itRule).left[0].id == nu.id) {
                rules.erase(itRule);
                break;
            }
        }
    }
}

void Grammar::recursiveInsertUnused(vector<Symbol> &unused, const Symbol& nt) {
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = rules[nt.id].right.begin();itRight != rules[nt.id].right.end()-1; itRight++) {
        unused.push_back((*itRight).first[1]);
        recursiveInsertUnused(unused, (*itRight).first[1]);
    }
}

void Grammar::normalizeProbs() {
    vector<Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        double freq = itRule->freq();
        for(itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++)
            itRight->second.first = itRight->second.first/freq;
    }

}

void Grammar::sampleRegularRules(
        vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> &vr,
        const std::vector<Symbol>& w) {
    Rule actualState = rules[0];
    pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> prod;
    for (const auto& s: w) {
        double total = 0.0;
        vector<double> probsTerminal;
        for (unsigned int i = nTerminals*s.id; i < nTerminals*s.id + nNonTerminals;i++) {
            probsTerminal.push_back(total + actualState.right[i].second.first);
            total += actualState.right[i].second.first;
        }
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        int i;
        for (i = 0; i < probsTerminal.size(); i++) {
            if (p * total < probsTerminal[i])
                break;
        }

        prod = make_pair(actualState.left, actualState.right[nTerminals*s.id+i]);
        vr.push_back(prod);
        actualState = rules[actualState.right[i].first[1].id];
    }
    prod = make_pair(actualState.left, actualState.right[actualState.right.size()-1]);
    vr.push_back(prod);
}

void Grammar::collapsedGibbsSamplePFA(int iterations) {
    //srand((unsigned) time(nullptr));
    int sentencesLength = words.size();
    for (int j = 0; j < iterations; j++) {
        parseTreesVec.clear();
        for (int i = 0; i< sentencesLength; i++) {
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
            sampleRegularRules(vProds,words[i]);
            std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTree;
            pairTree.first = words[i];
            pairTree.second = vProds;
            parseTreesVec.push_back(pairTree);
        }
        std::vector<Symbol> a;
        calculateNewThetaVec(a);
    }
}

void Grammar::generateNGramRules() {
    int lastSmallRuleI = 0;
    for (int i = 1; i < nNonTerminals; i++) {
        lastSmallRuleI += pow(terminals.size(), i);
    }

    for (int i = 0; i <= lastSmallRuleI; i++) {
        vector<Symbol> left;
        left.push_back(nonterminals[i]);
        vector<std::pair<std::vector<Symbol>, std::pair<double, double>>> right;
        for (const auto& t: terminals) {
            std::pair<std::vector<Symbol>, std::pair<double, double>> rightHS;
            rightHS.first.push_back(t);
            rightHS.first.push_back(nonterminals[nTerminals*i+t.id+1]);
            rightHS.second.first = 0;
            rightHS.second.second = ALFA;
            right.push_back(rightHS);
        }
        std::pair<std::vector<Symbol>, std::pair<double, double>> rightHS;
        //rightHS.second.first = 1.0/((nTerminals+1)*1.0);
        rightHS.second.first = 0;
        rightHS.second.second = ALFA;
        right.push_back(rightHS);
        Rule r = Rule(left, right);
        rules.push_back(r);
    }

    for (int i = lastSmallRuleI+1; i < nonterminals.size(); i++) {
        vector<Symbol> left;
        left.push_back(nonterminals[i]);
        vector<std::pair<std::vector<Symbol>, std::pair<double, double>>> right;
        for (const auto& t: terminals) {
            std::pair<std::vector<Symbol>, std::pair<double, double>> rightHS;
            rightHS.first.push_back(t);
            rightHS.first.push_back(nonterminals[lastSmallRuleI+1+(((i-lastSmallRuleI-1)% ((nonterminals.size()-lastSmallRuleI-1)/nTerminals))*nTerminals+t.id)]);
            //rightHS.second.first = 1.0/((nTerminals+1)*1.0);
            rightHS.second.first = 0;
            rightHS.second.second = ALFA;
            right.push_back(rightHS);
        }
        std::pair<std::vector<Symbol>, std::pair<double, double>> rightHS;
        rightHS.second.first = 0;
        rightHS.second.second = ALFA;
        right.push_back(rightHS);
        Rule r = Rule(left, right);
        rules.push_back(r);
    }

}


void Grammar::generateNGramNonTerminals() {
    int nIds = 0;
    Symbol nt = Symbol("NTI", nIds, false, false);
    nonterminals.push_back(nt);
    nt.name = "NT";
    for (int i = 1; i <= nNonTerminals; i++)
        recursiveAddTerminalToNT(nt, i, nIds);
}

void Grammar::recursiveAddTerminalToNT(Symbol nt, int n, int &nIds) {
    if (n ==0) {
        nonterminals.push_back(nt);
    } else {
        for (const auto& t: terminals) {
            nt.name += t.name;
            nt.id = ++nIds;
            recursiveAddTerminalToNT(nt, n-1,nIds);
            nt.name = nt.name.substr(0, nt.name.size()-1);
        }
    }
}

void Grammar::trainNGram() {
    for (auto w: words) {
        vector<Symbol> nGram;
        Symbol finalLHS = Symbol("NTI", 0, false, false);
        for(int i = 0; i < w.size(); i++) {
            if (i == nNonTerminals ) break;
            Symbol nt = Symbol("NTI", 0, false, false);
            if (i >=1)
                nt.name = "NT";
            finalLHS.name = "NT";
            for (const auto& s: nGram) {
                nt.name += s.name;
                finalLHS.name += s.name;
            }
            addNGramRuleFrequency(nt, w[i]);
            finalLHS.name += w[i].name;
            nGram.push_back(w[i]);
        }
        for(int i = nNonTerminals; i < w.size(); i++) {
            Symbol nt = Symbol("NT", 0, false, false);
            for (const auto& s: nGram)
                nt.name += s.name;
            addNGramRuleFrequency(nt, w[i]);
            finalLHS.name = nt.name;
            nGram.erase(nGram.begin());
            nGram.push_back(w[i]);
        }
        Symbol final = Symbol("", 0, true, false);
        addNGramRuleFrequency(finalLHS, final);
    }
    normalizeProbs();
}

void Grammar::addNGramRuleFrequency(const Symbol& lhs, const Symbol& nextSymbol) {
    vector<Rule>::iterator itRule;

    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        if (itRule->left[0].name == lhs.name) {
            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            bool flag = false;
            for (itRight = itRule->right.begin(); itRight != itRule->right.end()-1; itRight++){
                if (itRight->first[0].name == nextSymbol.name) {
                    itRight->second.first += 1.0;
                    flag = true;
                    break;
                }
            }
            if(flag) {
                break;
            } else {
                (itRule->right.end()-1)->second.first += 1.0;
                break;
            }
        }
    }
}

Grammar::Grammar(Grammar const &grammar) {

}




