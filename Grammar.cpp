//
// Created by henrique on 02/04/19.
//

#include <iostream>
#include <algorithm>
#include <tuple>
#include <cmath>
#include "Grammar.h"

Grammar::Grammar(const std::vector<Symbol> &terminals, const std::vector<Symbol> &nonterminals,
                 const std::vector<Rule> &rules, const Symbol &start) : terminals(terminals),
                                                                        nonterminals(nonterminals), rules(rules),
                                                                        start(start) {}

Grammar::Grammar(const std::vector<Symbol> &terminals, const std::vector<Symbol> &nonterminals,
                 const std::vector<Rule> &rules, const Symbol &start, int maxRightHandSize, int maxProductionRules,
                 int nTerminals, int nNonTerminals) : terminals(terminals), nonterminals(nonterminals), rules(rules),
                                                      start(start), maxRightHandSize(maxRightHandSize),
                                                      maxProductionRules(maxProductionRules), nTerminals(nTerminals),
                                                      nNonTerminals(nNonTerminals) {}

Grammar::Grammar(const std::vector<Symbol> &terminals, std::pair<int, int> contextSize,
                 int nNonTerminals, std::vector<std::vector<Symbol>> words, int trainingM) : terminals(terminals), contextSize(contextSize), nNonTerminals(nNonTerminals), words(words) {
    nTerminals = terminals.size();
    generateNonTermnals();
    generateRulesCNF();
    start = rules[0].left.front();


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
    std::string rulesStr = R"()";
    std::vector<Rule>::iterator itr;
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        rulesStr += R"(
                    )" +(*itr).ruleToStr();
    }
    return rulesStr;
}

void Grammar::generateNonTermnals() {

    for (int i = 0; i< nNonTerminals; i++) {
        nonterminals.push_back(Symbol("NT" + std::to_string(i), i, false, false));
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
    printGrammar();
}


double*** Grammar::CYKProb(std::string w) {
    //std::cout <<"IT for " << w << std::endl;
    double ***p = new double**[nonterminals.size()];
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
                if (itLeft->terminal == false && itLeft->context == false)
                    nonterminal = (*itLeft).clone();
                break;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal == true && itRightS->context == false) {
                        if(itRightS->name.compare(w.substr(i,1)) == 0) {
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
                        if (itLeft->terminal == false && itLeft->context == false)
                            nonterminal = (*itLeft).clone();
                        break;
                    }
                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (itRightS->terminal == false && itRightS->context == false) {
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

void Grammar::printInsideTable(double ***p, int wSize) {
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

void Grammar::sampleParseTree(std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>> &vr, Rule r, std::string w, double ***insideTable, int i, int k) {
    int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    if (jRange > 0) {
        //std::cout << w << std::endl;
        double *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).terminal == false && (*itSymbol).context == false) {
                    for (int l = 0; l < jRange; l++) {
                        int indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[(*itSymbol).id][i][i + l];
                        itSymbol++;
                        double pCjk = insideTable[(*itSymbol).id][i + l + 1][k];
                        int indexC = (*itSymbol).id * jRange;
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
        double p = (double) rand()/RAND_MAX;
        //std::cout << "Free: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        int j = 0;

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
            if ((*itRight).first[0].name.compare(s) == 0) {
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



void Grammar::calculateNewTheta(std::string w) {
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

int Grammar::calculateProductonCounts(std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>> production, std::string w) {
    int productionCount = 0;

    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol>,std::pair<std::vector<Symbol>,std::pair<double, double>>>>>>::iterator itParseTree;
    for (itParseTree = parseTrees.begin(); itParseTree != parseTrees.end(); itParseTree++) {
        if (w.compare((*itParseTree).first) != 0) {
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

std::vector<std::pair<double, int>> Grammar::calculateRuleFrequence(Rule r, std::string w) {
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

double Grammar::pTiTiMinus1 (std:: string w) {
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
    srand((unsigned) time(NULL));
    int sentencesLength = words.size();
    for (int i = 0; i< sentencesLength; i++) {
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
        double ***insideTable;
        insideTable = CYKProbVec(words[i]);
        sampleParseTreeVec(vProds, rules[0], words[i], insideTable, 0, words[i].size()-1);
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTree;
        pairTree.first = words[i];
        pairTree.second = vProds;
        parseTreesVec.push_back(pairTree);
    }
    int countAcceptedTress = 0;
    bool acceptTree = true;
    for (int j = 0; j < iterations; j++) {
        double ***insideTable;
        int i = rand()%sentencesLength;
        insideTable = CYKProbVec(words[i]);
        //printInsideTable(insideTable, sentences[i].size());
        if (j%100 == 0)
            std::cout << "iteration: " << j << " Tree " << i <<" Perplexity: "<< perplexity(words, true) << " Accepted Tress: " << countAcceptedTress << std::endl;
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
        double r = (double) rand()/RAND_MAX;

        if (!(funcA < r))
            parseTreesVec[i] = tpaux;
        else {
            acceptTree = true;
            countAcceptedTress++;
            std::cout << "Accepted tree in iteration: " << j << " Perplexity: " << perplexity(words, true) << std::endl;
            std::cout << " Pti'ti-1 = " << pTiLineTiMinis1 << " Ptiti-1 = "  << pTiTiMinis1<< " Ptiwi = " <<  pTiWi  << " Pti'wi = " << pTiLineWi << std::endl;
            std::cout << " ratio = " << (pTiLineTiMinis1*pTiWi)/(pTiTiMinis1*pTiLineWi) << ", funcA = " << funcA << ", r = " << r << std::endl;
        }
    }
    printGrammar();
}

void Grammar::metropolisHastingsPCSG(int iterations) {
    srand((unsigned) time(NULL));
    int sentencesLength = words.size();
    for (int i = 0; i< sentencesLength; i++) {
        actualProduction.clear();
        actualProduction.push_back(nonterminals[0]);
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
        double ****insideTableKL;
        insideTableKL = CYKProbKLVec(words[i]);
        sampleParseTreeKLVec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size()-1);
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTreeVec;
        pairTreeVec.first = words[i];
        pairTreeVec.second = vProds;
        parseTreesVec.push_back(pairTreeVec);
    }
    for (int j = 0; j < iterations; j++) {
        int i = (int) ((rand()/RAND_MAX) * (sentencesLength-1));
        double ****insideTableKL;
        insideTableKL = CYKProbKLVec(words[i]);
        //printInsideTableKL(insideTableKL,sentences[i].size());
        actualProduction.clear();
        actualProduction.push_back(nonterminals[0]);
        if (j%1 == 0)
            std::cout << "iteration: " << j << " Tree " << i << " Perplexity: "<< perplexityKL(words, true) <<std::endl;
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;

        //std::cout << "SENSITIVE" << std::endl;
        calculateNewThetaVec(words[i]);
        updateParseTressTheta();

        insideTableKL = CYKProbKLVec(words[i]);
        //printInsideTableKL(insideTableKL,sentences[i].size());

        vProds = parseTreesVec[i].second;
        double pTiWi = probTree(vProds)/insideTableKL[0][0][words[i].size()-1][0] ;
        double pTiTiMinis1 = pTiTiMinus1Vec(words[i]);
        vProds.clear();

        sampleParseTreeKLVec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size()-1);
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> tp;
        std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> tpaux;
        tpaux = parseTreesVec[i];
        tp.first = parseTreesVec[i].first,
        tp.second = vProds;
        parseTreesVec[i] = tp;
        double pTiLineWi =  probTree(vProds) / insideTableKL[0][0][words[i].size()-1][0];
        double pTiLineTiMinis1 = pTiTiMinus1Vec(words[i]);
        double funcA = std::min(1.0,(pTiLineTiMinis1*pTiWi)/(pTiTiMinis1*pTiLineWi));
        if (!(funcA < (double) rand()/RAND_MAX)) {
            parseTreesVec[i] = tpaux;
        }
        else
            std::cout << "Accepted tree in iteration: " << j << " Perplexity: " << perplexityKL(words, true) << std::endl;
    }
    printGrammar();
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
                        r.generatePiorDirichlet(ALFA,0);
                        r.updateProbDirichletTheta();
                        rules.push_back(r);
                    }
                }
            }
        }
    }


}

double ****Grammar::CYKProbKL(std::string w) {
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
    double ****p = new double***[nonterminals.size()*sizeLeftContext];
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
                if (itLeft->terminal == false && itLeft->context == false) {
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
                    if (itRightS->terminal == true && itRightS->context == false) {
                        if(itRightS->name.compare(w.substr(i,1)) == 0) {
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
                        if (itLeft->terminal == false && itLeft->context == false) {
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
                            if (itRightS->terminal == false && itRightS->context == false) {
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

int Grammar::convertContextToID(int side, std::vector<Symbol> context) {
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
        id += (*itSymbol).id * pow(nSymbols, position);
        position++;
    }
    return id;
}

void Grammar::sampleParseTreeKL(
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> &vr,
        Rule r, std::string w, double ****insideTable, int i, int k) {
    int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    Symbol nonterminal =  Symbol("", 0,false);
    std::vector<Symbol>::iterator itLeft;
    std::vector<Symbol> leftContext;
    for (itLeft = r.left.begin(); itLeft != r.left.end(); itLeft++) {
        if (itLeft->terminal == false && itLeft->context == false) {
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
        double *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).terminal == false && (*itSymbol).context == false) {
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
                        int indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                                [i][i + l][convertContextToID(1,rightContextLine)];
                        itSymbol++;
                        double pCjk = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                                [i + l + 1][k][convertContextToID(1,rightContext)];
                        int indexC = (*itSymbol).id * jRange;
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
        double p = (double) rand()/RAND_MAX;
        //std::cout << "Sensitive: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        int j = 0;


        Symbol rightB = Symbol("", 0,false);;
        Symbol rightC = Symbol("", 0,false);;

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

        int lc , rc;
        lc = ((int) rand()%(leftContext2.size()+1));
        rc = (int) rand()%(rightContext2.size()+1);
        for (int i = 0; i < lc; i++)
            leftContext2.erase(leftContext2.begin());
        for (int i = 0; i < rc; i++)
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
        lc = (int) rand()%(leftContext2.size()+1);
        rc = (int) rand()%(rightContext2.size()+1);
        for (int i = 0; i < lc; i++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i = 0; i < rc; i++)
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
        int lc , rc;
        lc = (int) rand()%(leftContext2.size()+1);
        rc = (int) rand()%(rightContext2.size()+1);
        for (int i = 0; i < lc; i++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i = 0; i < rc; i++)
            rightContext2.pop_back();

        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w.substr(i,jRange+1);
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            std::vector<Symbol>::iterator itSymbol;
            bool flagContext = false;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).name.compare(s) == 0 && !(*itSymbol).context ) {
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
                        if (itRight->terminal == false && itRight->context == false) {
                            newProd.push_back(*itRight);
                            itRight++;
                            newProd.push_back(*(itRight));
                            appliedProd = true;
                        }
                        if (itRight->terminal == true && itRight->context == false) {
                            newProd.push_back(*itRight);
                            appliedProd = true;
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
        std::vector<Symbol>::iterator itRuleLeft = (*itRule).left.begin();
        bool flag = true;
        for (itLHS = lhs.begin(); itLHS != lhs.end(); itLHS++) {

            if ((*itRuleLeft).name.compare((*itLHS).name) != 0)
                flag = false;
            if (itRuleLeft+1 != (*itRule).left.end())
                itRuleLeft++;
        }
        if (flag)
            return (*itRule);
    }
}

void Grammar::getActualContext(std::vector<Symbol> &leftContext, std::vector<Symbol> &rightContext) {
    leftContext.clear();
    rightContext.clear();
    std::vector<Symbol>::iterator itLeft;
    Symbol aux = Symbol("",0,true);
    for (itLeft = actualProduction.begin(); itLeft != actualProduction.end(); itLeft++) {
        if (itLeft->terminal == false && itLeft->context == false) {
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
        Rule r, std::vector<Symbol> w, double ***insideTable, int i, unsigned long k) {
    int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    if (jRange > 0) {
        //std::cout << w << std::endl;
        double *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).terminal == false && (*itSymbol).context == false) {
                    for (int l = 0; l < jRange; l++) {
                        int indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[(*itSymbol).id][i][i + l];
                        itSymbol++;
                        double pCjk = insideTable[(*itSymbol).id][i + l + 1][k];
                        int indexC = (*itSymbol).id * jRange;
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
        double p = (double) rand()/RAND_MAX;
        //std::cout << "Free: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        int j = 0;

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
            if ((*itRight).first[0].name.compare(s) == 0) {
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
    double ***p = new double**[nonterminals.size()];
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
                if (itLeft->terminal == false && itLeft->context == false)
                    nonterminal = (*itLeft).clone();
                break;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal == true && itRightS->context == false) {
                        if(itRightS->name.compare(w[i].name) == 0) {
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
                        if (itLeft->terminal == false && itLeft->context == false)
                            nonterminal = (*itLeft).clone();
                        break;
                    }
                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (itRightS->terminal == false && itRightS->context == false) {
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

void Grammar::calculateNewThetaVec(std::vector<Symbol> w) {
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
        std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>> production,
        std::vector<Symbol> w) {
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
    std::vector<Symbol>::iterator itW2 = w2.begin();
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

double Grammar::pTiTiMinus1Vec(std::vector<Symbol> w) {
    std::vector<Rule>::iterator itRule;
    double result = 1;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<Symbol> a;
        result *=  cConstant(calculateRuleFrequenceVec((*itRule),  a))/cConstant(calculateRuleFrequenceVec((*itRule), w));
        //std::cout << "   PTia = " << cConstant(calculateRuleFrequenceVec((*itRule),  a)) << " Pt-ia = " << cConstant(calculateRuleFrequenceVec((*itRule), w)) << std::endl;
        //Talvez essa passagem de parametro esteja errada

    }
}
std::vector<std::pair<double, int>> Grammar::calculateRuleFrequenceVec(Rule &r, std::vector<Symbol> w) {
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
    double ****p = new double***[nonterminals.size()*sizeLeftContext];
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
                if (itLeft->terminal == false && itLeft->context == false) {
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
                    if (itRightS->terminal == true && itRightS->context == false) {
                        if(itRightS->name.compare(w[i].name) == 0) {
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
                        if (itLeft->terminal == false && itLeft->context == false) {
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
                            if (itRightS->terminal == false && itRightS->context == false) {
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

void Grammar::sampleParseTreeKLVec(
        std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>, std::pair<double, double>>>> &vr,
        Rule r, std::vector<Symbol> w, double ****insideTable, int i, int k) {
    int jRange = k - i;
    int rightRange = nNonTerminals*nNonTerminals;

    Symbol nonterminal =  Symbol("", 0,false);
    std::vector<Symbol>::iterator itLeft;
    std::vector<Symbol> leftContext;
    for (itLeft = r.left.begin(); itLeft != r.left.end(); itLeft++) {
        if (itLeft->terminal == false && itLeft->context == false) {
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
        double *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).terminal == false && (*itSymbol).context == false) {
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
                        int indexB = (*itSymbol).id * nNonTerminals * jRange;
                        double pBij = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                        [i][i + l][convertContextToID(1,rightContextLine)];
                        itSymbol++;
                        double pCjk = insideTable[contextAmount.first*(*itSymbol).id+ convertContextToID(0,leftContext)]
                        [i + l + 1][k][convertContextToID(1,rightContext)];
                        int indexC = (*itSymbol).id * jRange;
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
        double p = (double) rand()/RAND_MAX;
        //std::cout << "Sensitive: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        int j = 0;


        Symbol rightB = Symbol("", 0,false);;
        Symbol rightC = Symbol("", 0,false);;

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

        int lc , rc;
        lc = ((int) rand()%(leftContext2.size()+1));
        rc = (int) rand()%(rightContext2.size()+1);
        for (int i = 0; i < lc; i++)
            leftContext2.erase(leftContext2.begin());
        for (int i = 0; i < rc; i++)
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
        lc = (int) rand()%(leftContext2.size()+1);
        rc = (int) rand()%(rightContext2.size()+1);
        for (int i = 0; i < lc; i++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i = 0; i < rc; i++)
            rightContext2.pop_back();

        lhs.clear();
        lhs.insert(lhs.begin(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightC);
        lhs.insert(lhs.begin(), rightContext2.begin(), rightContext2.end());
        sampleParseTreeKLVec(vr, findRuleByLHS(lhs), w, insideTable, i+j+1, k);

    } else {
        std::vector<Symbol> leftContext2;
        std::vector<Symbol> rightContext2;
        getActualContext(leftContext2, rightContext2);

        while (leftContext2.size() > contextSize.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > contextSize.second)
            rightContext2.pop_back();
        int lc , rc;
        lc = (int) rand()%(leftContext2.size()+1);
        rc = (int) rand()%(rightContext2.size()+1);
        for (int i = 0; i < lc; i++)
            leftContext2.erase(leftContext2.begin() +i);
        for (int i = 0; i < rc; i++)
            rightContext2.pop_back();
        std::pair<std::vector<Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w[i].name;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            std::vector<Symbol>::iterator itSymbol;
            bool flagContext = false;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).name.compare(s) == 0 && !(*itSymbol).context ) {
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

double Grammar::perplexity(std::vector<std::vector<Symbol>> testData, bool normalized) {
    double pXTheta = 1.0;
    for (auto w: testData) {
        double ***iTable  = CYKProbVec(w);
        double pX = iTable[0][0][w.size()-1] ;
        if (normalized) {
            double ***pN = CYKProbNVec(w);
            pX /= pN[0][0][w.size()-1];
        }
        pXTheta *= pX ;
    }
    return exp(-((1.0/testData.size()) * log(pXTheta)));
}

double Grammar::perplexityKL(std::vector<std::vector<Symbol>> testData, bool normalized) {
    double pXTheta = 1.0;
    for (auto w: testData) {
        double ****iTable  = CYKProbKLVec(w);
        double pX = iTable[0][0][w.size()-1][0] ;
        if (normalized) {
            double ****pN = CYKProbKLNVec(w);
            pX /= pN[0][0][w.size()-1][0];
        }
        pXTheta *= pX;
    }
    return exp(-((1.0/testData.size()) * log(pXTheta)));
}

double ***Grammar::CYKProbNVec(std::vector<Symbol> w) {
    double ***p = new double**[nonterminals.size()];
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
                if (itLeft->terminal == false && itLeft->context == false)
                    nonterminal = (*itLeft).clone();
                break;
            }

            std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal == true && itRightS->context == false) {
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
                        if (itLeft->terminal == false && itLeft->context == false)
                            nonterminal = (*itLeft).clone();
                        break;
                    }
                    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (itRightS->terminal == false && itRightS->context == false) {
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

double ****Grammar::CYKProbKLNVec(std::vector<Symbol> w) {
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
    double ****p = new double***[nonterminals.size()*sizeLeftContext];
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
                if (itLeft->terminal == false && itLeft->context == false) {
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
                    if (itRightS->terminal == true && itRightS->context == false) {
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
                        if (itLeft->terminal == false && itLeft->context == false) {
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
                            if (itRightS->terminal == false && itRightS->context == false) {
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

Grammar::~Grammar() {

}

void Grammar::gibbsSamplingPCFG(int iterations) {
    srand((unsigned) time(NULL));
    int sentencesLength = words.size();
    for (int j = 0; j < iterations; j++) {
        parseTreesVec.clear();
        for (int i = 0; i< sentencesLength; i++) {
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
            double ***insideTable;
            insideTable = CYKProbVec(words[i]);
            sampleParseTreeVec(vProds, rules[0], words[i], insideTable, 0, words[i].size()-1);
            std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTree;
            pairTree.first = words[i];
            pairTree.second = vProds;
            parseTreesVec.push_back(pairTree);
        }
        std::vector<Symbol> a;
        calculateNewThetaVec(a);
        updateParseTressTheta();
        if (j%5 == 0)
            std::cout << "iteration: " << j <<" Perplexity: "<< perplexity(words, true) << std::endl;
    }
}

void Grammar::gibbsSamplingPCSG(int iterations) {
    srand((unsigned) time(NULL));
    int sentencesLength = words.size();
    for (int j = 0; j < iterations; j++) {
        parseTreesVec.clear();
        for (int i = 0; i< sentencesLength; i++) {
            actualProduction.clear();
            actualProduction.push_back(nonterminals[0]);
            std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>> vProds;
            double ****insideTableKL;
            insideTableKL = CYKProbKLVec(words[i]);
            sampleParseTreeKLVec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size()-1);
            std::pair<std::vector<Symbol>, std::vector<std::pair<std::vector<Symbol>, std::pair<std::vector<Symbol>,std::pair<double, double>>>>> pairTreeVec;
            pairTreeVec.first = words[i];
            pairTreeVec.second = vProds;
            parseTreesVec.push_back(pairTreeVec);
        }
        std::vector<Symbol> a;
        calculateNewThetaVec(a);
        updateParseTressTheta();
        if (j%5 == 0)
            std::cout << "iteration: " << j <<" Perplexity: "<< perplexity(words, true) << std::endl;
    }
}
