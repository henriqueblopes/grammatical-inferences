//
// Created by henrique on 02/04/19.
//

#include "Grammar.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <langinfo.h>
#include <queue>
#include <sstream>
#include <stack>
#include <tuple>
#include <utility>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
using namespace std;


Grammar::Grammar::Grammar(const std::vector<Symbol::Symbol> &terminals, int nNonTerminals, std::vector<std::vector<Symbol::Symbol>> words,
                 enum grammar_type g_tp, pair<int,int> contextSize) : terminals(terminals), n_non_terminals(nNonTerminals), words(std::move(std::move(words))), g_tp(g_tp) {
    n_terminals = terminals.size();
    this->context_size = contextSize;
    if (g_tp == n_gram ) {
        contextSize = make_pair(0,0);
        generate_n_gram_non_terminals();
        generate_n_gram_rules();
    } else if (g_tp == pfa || g_tp == pcfg){
        contextSize = make_pair(0,0);
        generate_non_termnals();
        if(g_tp == pfa) {
            contextSize = make_pair(0,0);
            generate_rules_regular();
            start = rules[0].left.front();
        }
        else if (g_tp == pcfg) {
            this->context_size = make_pair(0, 0);
            generate_rules_cnf();
            start = rules[0].left.front();
        }
    } else {
        n_terminals = terminals.size();
        generate_non_termnals();
        generate_rules_cnf();
        start = rules[0].left.front();
    }

}


void Grammar::Grammar::print_grammar() {
    //std::cout << "Terminals: " << std::endl;
    std::vector<Symbol::Symbol>::iterator it;
    /*for(it = terminals.begin(); it != terminals.end(); it++) {
        std::cout << "\t";
        std::cout <<  (*it).name << std::endl;
    }

    std::cout << std::endl;*/

    /*std::cout << "NonTerminals: " << std::endl;
    for(it = non_terminals.begin(); it != non_terminals.end(); it++) {
        std::cout << "\t";
        std::cout <<  (*it).name <<std::endl;
    }

    std::cout << std::endl;*/

    std::vector<Rule::Rule>::iterator itr;
    std::cout << "Rules: " << std::endl;
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        std::cout << "\t";
        (*itr).print_rule();
    }

    std::cout << std::endl;
}

[[maybe_unused]] std::string Grammar::Grammar::grammar_to_str() {
    std::string g;
    g.append("Terminals: ");
    g.append("\n");
    std::vector<Symbol::Symbol>::iterator it;
    for(it = terminals.begin(); it != terminals.end(); it++) {
        g.append("\t");
        g.append((*it).name);
        g.append("\n");
    }

    g.append("\n");

    g.append("NonTerminals: ");
    g.append("\n");
    for(it = non_terminals.begin(); it != non_terminals.end(); it++) {
        g.append("\t");
        g.append((*it).name);
        g.append("\n");
    }
    g.append("\n");

    std::vector<Rule::Rule>::iterator itr;
    g.append("Rules: ");
    g.append("\n");
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        g.append("\t");
        g.append((*itr).rule_to_str_lalr());
    }
    g.append("\n");
    return g;
}

void Grammar::Grammar::print_rules() {
    std::vector<Rule::Rule>::iterator itr;
    std::cout << "Rules: " << std::endl;
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        std::cout << "\t";
        (*itr).print_rule();
    }

    std::cout << std::endl;
}

[[maybe_unused]] std::string Grammar::Grammar::rules_to_string () {
    std::string rulesStr;
    std::vector<Rule::Rule>::iterator itr;
    for(itr = rules.begin(); itr != rules.end(); itr++) {
        rulesStr += R"(
                    )" + (*itr).rule_to_str();
    }
    return rulesStr;
}

void Grammar::Grammar::generate_non_termnals() {

    for (size_t i = 0; i < n_non_terminals; i++) {
        non_terminals.emplace_back("NT" + std::to_string(i), i, false, false);
    }
}



void Grammar::Grammar::generate_permutation(std::vector<std::vector<Symbol::Symbol>> & permutations, std::vector<Symbol::Symbol> symbols, size_t size, std::vector<Symbol::Symbol> word, bool context) {
    if (size == 0 )
        permutations.push_back(word);
    else {

        std::vector<Symbol::Symbol>::iterator itSymbol;
        for (itSymbol = symbols.begin(); itSymbol != symbols.end(); itSymbol++) {
            Symbol::Symbol aux = (*itSymbol).clone();
            aux.context = context;
            word.push_back(aux);
            generate_permutation(permutations, symbols, size - 1, word, context);
            word.pop_back();
        }
    }
}

void Grammar::Grammar::train(training_method algorithm, int iterations) {
    switch (g_tp) {
        case n_gram:
            train_n_gram();
            break;
        case pfa:
            switch(algorithm) {
                case pfa_baum_welch:
                    baum_welch(iterations);
                    break;
                case pfa_collapsed_gibbs_sample:
                    collapsed_gibbs_sample_pfa(iterations);
                    break;
                case pfa_alergia:
                    alergia(0.01);
                    break;
                default:
                    baum_welch(iterations);
                    break;
            }
            break;
        case pcfg:
            switch(algorithm) {
                case pcfg_inside_outside:
                    inside_outside(iterations);
                    break;
                case pcfg_metropolis_hastings:
                    metropolis_hastings_pcfg(iterations);
                    break;
                case pcfg_gibbs_sampling:
                    gibbs_sampling_pcfg(iterations);
                    break;
                default:
                    inside_outside(iterations);
                    break;
            }
            break;
        case pcsg:
            switch(algorithm) {
                case pcsg_metropolis_hastings:
                    metropolis_hastings_pcsg(iterations);
                    break;
                case pcsg_gibbs_sampling:
                    gibbs_sampling_pcsg(iterations);
                    break;
                default:
                    metropolis_hastings_pcsg(iterations);
                    break;
            }
            break;
    }/*
    if (algorithm == 1) {
        context_size.first = context_size.second = 0;
        metropolis_hastings_pcfg(iterations);
    }
    else if (algorithm == 3) {
        context_size.first = context_size.second = 0;
        gibbs_sampling_pcfg(iterations);
    }
    else if (algorithm == 4) {
        gibbs_sampling_pcsg(iterations);
    }
    else
        metropolis_hastings_pcsg(iterations);*/
    //print_grammar();
}



void Grammar::Grammar::print_inside_table(double ***p, int wSize) const {
    cout << "IutsideTable" << endl;
    for (unsigned long i = 0; i < non_terminals.size(); i++) {
        std::cout <<"Nonterminal: " << i << "\n";
        for (int j = 0; j < wSize; j++) {
            for (int k = 0; k < wSize; k++) {
                std::cout << p[i][j][k] << " ";
            }
            std::cout << "\n";
        }
    }
}

void Grammar::Grammar::print_outside_table(const std::vector<std::vector<std::vector<double>>>& p) {
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

/*
void Grammar::printInsideTableKL(double ****p, int wSize) {
    int sizeLeftContext = 0;
    int sizeRightContext = 0;
    for (unsigned long i = 0; i <= context_size.first; i++) {
        std::vector<Symbol::Symbol> word;
        std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();
    }
    for (unsigned long i = 0; i <= context_size.second; i++) {
        std::vector<Symbol::Symbol> word;
        std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
        generate_permutation(permutationsTerminals, non_terminals, i, word, true);
        sizeRightContext += permutationsTerminals.size();
    }
    for (unsigned long i = 0; i < non_terminals.size()*sizeLeftContext; i++) {
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
*/

void Grammar::Grammar::sample_parse_tree(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> &vr, Rule::Rule r, const std::string& w, double ***inside_table, size_t i, size_t k) {
    size_t jRange = k - i;
    size_t rightRange = n_non_terminals * n_non_terminals;

    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    for (unsigned long l = 0; l < jRange; l++) {
                        size_t indexB;
                        indexB = (*itSymbol).id * n_non_terminals * jRange;
                        double pBij = inside_table[(*itSymbol).id][i][i + l];
                        itSymbol++;
                        double pCjk = inside_table[(*itSymbol).id][i + l + 1][k];
                        size_t indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) / inside_table[(r.left[0].id)][i][k];
                        itSymbol--;
                        /*std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first << " * " <<
                        pBij << " ("<< "IT["<<(r.left[0].id)<<"]["<<i<<"]["<<i+l<<"])"<< " * " <<
                        pCjk << " ("<< "IT["<<(r.left[0].id)<<"]["<<i+l+1<<"]["<<k<<"])" << " / " << inside_table[(r.left[0].id)][i][k]
                        <<"IT["<<(r.left[0].id)<<"]["<<i<<"]["<<k<<"]"<<std::endl;*/
                    }
                    break;
                }

            }
        }
        size_t l = 0;
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
        size_t j = 0;
        delete[]jbc;
        //std::vector<Symbol::Symbol> rightProduction;
        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    if ((*itSymbol).id == l / (n_non_terminals * jRange)) {
                        itSymbol++;
                        if(itSymbol == (*itRight).first.end())
                            break;
                        if((*itSymbol).id == (l%(n_non_terminals * jRange)) / jRange) {
                            itSymbol--;
                            rightProduction.first.push_back((*itSymbol));
                            j = l - n_non_terminals * jRange * (*itSymbol).id;
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
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        sample_parse_tree(vr, rules[production.second.first[0].id], w, inside_table, i, i + j);
        sample_parse_tree(vr, rules[production.second.first[1].id], w, inside_table, i + j + 1, k);

    } else {
        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w.substr(i,jRange+1);
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            //TO DO:: Dar um jeito de não hardecodá o left[0]
            if ((*itRight).first[0].name == s) {
                rightProduction.first.push_back((*itRight).first[0]);
                rightProduction.second.first = (*itRight).second.first;
                rightProduction.second.second = (*itRight).second.second;
            }
        }
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        return;
    }

}
void Grammar::Grammar::calculate_new_theta_vec_opt(int i) {
    p_ti_minus_1_frequence(i);
    std::vector<Rule::Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        double nonterminalTotal = 0;
        std::vector<std::pair<double, int>>::iterator iTrFrequence;
        for (iTrFrequence = (*itRule).rule_frequence.begin(); iTrFrequence != (*itRule).rule_frequence.end(); iTrFrequence++) {
            nonterminalTotal += iTrFrequence->second + iTrFrequence->first;
        }
        iTrFrequence = itRule->rule_frequence.begin();
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            (*itRight).second.first = ((*iTrFrequence).second + (*itRight).second.second )/nonterminalTotal;
            iTrFrequence++;
        }
    }
    p_ti_minus_1_plus_frequence(i);
}


int Grammar::Grammar::calculate_producton_counts(const std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>& production, const std::string& w) {
    int productionCount = 0;

    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>>>::iterator itParseTree;
    for (itParseTree = parse_trees.begin(); itParseTree != parse_trees.end(); itParseTree++) {
        if (w != (*itParseTree).first) {
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>::iterator itProduction;
            for (itProduction = (*itParseTree).second.begin(); itProduction != (*itParseTree).second.end(); itProduction++) {
                if(equal_productions((*itProduction), production)) {
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

std::vector<std::pair<double, int>> Grammar::Grammar::calculate_rule_frequence(Rule::Rule r, const std::string& w) {
    std::vector<std::pair<double, int>> ruleFrequence;
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production;
        production.first = r.left;
        production.second.first = (*itRight).first;
        std::pair<double, int> productionFrequence;
        productionFrequence.first = (*itRight).second.second;
        productionFrequence.second = calculate_producton_counts(production, w);
        ruleFrequence.push_back(productionFrequence);
    }
    return ruleFrequence;
}

double Grammar::Grammar::c_constant(std::vector<std::pair<double, int>> ruleFrequence) {
    double numerator = 1.0;
    double denominator = 0.0;
    std::vector<std::pair<double, int>>::iterator itFrequence;
    for (itFrequence = ruleFrequence.begin(); itFrequence != ruleFrequence.end(); itFrequence++) {
        numerator = numerator * tgamma((*itFrequence).first + (*itFrequence).second);
        denominator = denominator + (*itFrequence).first + (*itFrequence).second;
    }
    return numerator/tgamma(denominator);
}


double Grammar::Grammar::prob_tree(std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> tree) {
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>::iterator itTree;
    double prob = 1.0;
    for (itTree = tree.begin(); itTree != tree.end(); itTree++) {
        prob *= (*itTree).second.second.first;
    }
    return prob;
}

bool Grammar::Grammar::equal_productions (std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> prod_a, std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> prod_b) {
    std::vector<Symbol::Symbol>::iterator itA;
    std::vector<Symbol::Symbol>::iterator itB;
    bool flag = true;
    itB = prod_b.first.begin();
    for (itA = prod_a.first.begin(); itA != prod_a.first.end(); itA++) {
        if(!(*itA).equal_symbol((*itB))) {
            flag = false;
            break;
        }
        itB++;
        if(itB == prod_b.first.end() && (itA + 1) != prod_a.first.end()) {
            flag = false;
            break;
        }

    }
    if(itB != prod_b.first.end())
        flag = false;
    if(flag) {
        itB = prod_b.second.first.begin();
        for (itA = prod_a.second.first.begin(); itA != prod_a.second.first.end(); itA++) {
            //retirar a condição abaixo caso dê merda em outros métodos
            if (itB == prod_b.second.first.end()) {
                flag = false;
                break;
            }
            if(!(*itA).equal_symbol((*itB))) {
                flag = false;
                break;
            }
            itB++;
            if(itB == prod_b.second.first.end() && (itA + 1) != prod_a.second.first.end()) {
                flag = false;
                break;
            }
        }
        if(itB != prod_b.second.first.end())
            flag = false;
    }
    return flag;
}

void Grammar::Grammar::update_parse_tress_theta () {
    std::vector<std::pair<std::string, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>>>::iterator itParseTrees;

    for (itParseTrees = parse_trees.begin(); itParseTrees != parse_trees.end(); itParseTrees++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>::iterator itProduction;
        for (itProduction = (*itParseTrees).second.begin(); itProduction != (*itParseTrees).second.end(); itProduction++) {
            std::vector<Rule::Rule>::iterator itRule;
            bool flagProduction = false;
            for(itRule = rules.begin(); itRule != rules.end(); itRule++) {
                std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
                for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                    std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production;
                    production.first = (*itRule).left;
                    production.second = (*itRight);
                    if (equal_productions((*itProduction), production)) {
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

/*
void Grammar::printProduction(std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> tree) {
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>::iterator  itProduction;
    std::cout << "\nProductions:";
    for (itProduction = tree.begin(); itProduction != tree.end(); itProduction++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itProduction).first.begin(); itSymbol != (*itProduction).first.end(); itSymbol++)
                std::cout << (*itSymbol).name << " ";
            std::cout << " -> ";
            for (itSymbol = (*itProduction).second.first.begin(); itSymbol != (*itProduction).second.first.end(); itSymbol++)
                std::cout << (*itSymbol).name << " ";
            std::cout << "\n";
    }
}

void Grammar::printOneProduction(std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> prd) {
    std::vector<Symbol::Symbol>::iterator itSymbol;
    for (itSymbol = prd.first.begin(); itSymbol != prd.first.end(); itSymbol++)
        std::cout << (*itSymbol).name << " ";
    std::cout << " -> ";
    for (itSymbol = (prd).second.first.begin(); itSymbol != (prd).second.first.end(); itSymbol++)
        std::cout << (*itSymbol).name << " ";
    std::cout << "\n";
}
*/


void Grammar::Grammar::metropolis_hastings_pcfg(int iterations) {
    //srand((unsigned) time(nullptr));
    size_t sentencesLength = words.size();
    for (size_t i = 0; i< sentencesLength; i++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> vProds;
        double ***insideTable;
        insideTable = cyk_prob_vec(words[i]);
        sample_parse_tree_vec(vProds, rules[0], words[i], insideTable, 0, words[i].size() - 1);
        free_inside_table(insideTable, words[i].size());
        std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> pairTree;
        pairTree.first = words[i];
        pairTree.second = vProds;
        parse_trees_vec.push_back(pairTree);
    }
    int countAcceptedTress = 0;
    //bool acceptTree = true;
    for (int j = 0; j < iterations; j++) {
        double ***insideTable;

        std::mt19937 mt(rd());
        std::uniform_int_distribution<size_t> dist(0, sentencesLength-1);
        size_t p = dist(mt);
        size_t i = p;
        //int i = rand()%sentencesLength;

        //print_inside_table(insideTable, sentences[i].size());
        if (j%(iterations/10) == 0) {
            //std::pair<double,  double> pMpP = perplexity(words, true);
            std::cout << "iteration: " << j << " Tree " << i <<" Accepted Tress: " << countAcceptedTress << std::endl;
        }

        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> vProds;
        //std::cout << "FREE" << std::endl;
        //if (acceptTree) {
        calculate_new_theta_vec(words[i]);
        update_parse_tress_theta();
            //acceptTree = false;
        //}
        insideTable = cyk_prob_vec(words[i]);
        //print_inside_table(insideTable, sentences[i].size());
        vProds = parse_trees_vec[i].second;
        double pTiWi = prob_tree(vProds) / insideTable[0][0][words[i].size() - 1] ;
        double pTiTiMinis1 = p_ti_ti_minus_1_vec(words[i]);
        vProds.clear();

        sample_parse_tree_vec(vProds, rules[0], words[i], insideTable, 0, words[i].size() - 1);
        std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> tp;
        std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> tpaux;
        tpaux = parse_trees_vec[i];
        tp.first = parse_trees_vec[i].first,
                tp.second = vProds;
        parse_trees_vec[i] = tp;
        double pTiLineWi = prob_tree(vProds) / insideTable[0][0][words[i].size() - 1];
        double pTiLineTiMinis1 = p_ti_ti_minus_1_vec(words[i]);
        double funcA = std::min(1.0,(pTiLineTiMinis1*pTiWi)/(pTiTiMinis1*pTiLineWi));
        std::mt19937 md(rd());
        std::uniform_real_distribution<double> dd(0, 1);
        double r = dd(md);
        free_inside_table(insideTable, words[i].size());

        if (funcA >= r)
            parse_trees_vec[i] = tpaux;
        else {
//            acceptTree = true;
            countAcceptedTress++;
            /*std::pair<double,  double> pMpP = perplexity(words, true);
            std::cout << "  Accepted tree in iteration: " << j << " - PerplexityM: " << pMpP.first << " = PerplexityP: " << pMpP.second<< std::endl;
            std::cout << "    Pti'ti-1 = " << pTiLineTiMinis1 << " Ptiti-1 = "  << pTiTiMinis1<< " Ptiwi = " <<  pTiWi  << " Pti'wi = " << pTiLineWi << std::endl;
            std::cout << "    ratio = " << (pTiLineTiMinis1*pTiWi)/(pTiTiMinis1*pTiLineWi) << ", funcA = " << funcA << ", r = " << r << std::endl;*/
        }
    }
    //print_grammar();
}

void Grammar::Grammar::metropolis_hastings_pcsg(int iterations) {
    std::chrono::duration<double> totalTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> insideTTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> insideTTimeMet = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> iterationTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> newThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> updateThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> perplexityTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> remainingTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    //srand((unsigned) time(nullptr));
    size_t sentencesLength = words.size();
    for (size_t i = 0; i< sentencesLength; i++) {
        if (sentencesLength >= 10)
            if (i%(sentencesLength/10) ==0)
                std::cout << 100*(i/(1.0*sentencesLength))<< "% of trees parsed" << std::endl;
        auto startIt = std::chrono::system_clock::now();
        actual_production.clear();
        actual_production.push_back(non_terminals[0]);
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> vProds;
        double ****insideTableKL ;
        auto startInsideTable = std::chrono::system_clock::now();
        insideTableKL = cyk_prob_kl_vec(words[i]);
        insideTTime += std::chrono::system_clock::now() - startInsideTable;
        sample_parse_tree_kl_vec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size() - 1);
        free_inside_table_kl(insideTableKL, words[i].size());
        std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> pairTreeVec;
        pairTreeVec.first = words[i];
        pairTreeVec.second = vProds;
        parse_trees_vec.push_back(pairTreeVec);
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
        //double ****iTableN;
        actual_production.clear();
        actual_production.push_back(non_terminals[0]);
        if (iterations >= 5)
            if (j%(iterations/5) == 0) {
            //auto startPerplexityTime = std::chrono::system_clock::now();
            //std::pair<double,  double> pMpP = perplexity_kl(words, true);
            //perplexityTime += std::chrono::system_clock::now() - startPerplexityTime;
            //std::cout << "   iteration: " << j << " Tree " << i <<" - Perplexity: "<< pMpP.first << " - PerplexityN: "<< pMpP.second <<" Accepted Tress: " << countAcceptedTress << " PTime: "<< perplexityTime.count() <<std::endl;
                std::cout << "   iteration: " << j << " Tree " << i << " Accepted Tress: " << countAcceptedTress << std::endl;
        }
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> vProds;

        //std::cout << "SENSITIVE" << std::endl;
        auto startNewTheta = std::chrono::system_clock::now();
        calculate_new_theta_vec_opt(i);
        newThetaTime += std::chrono::system_clock::now() - startNewTheta;

        update_parse_tress_theta();

        auto startInsideTable = std::chrono::system_clock::now();
        insideTableKL = cyk_prob_kl_vec(words[i]);
        insideTTimeMet += std::chrono::system_clock::now() - startInsideTable;

        auto startRemaining = std::chrono::system_clock::now();
        vProds = parse_trees_vec[i].second;

        double pTiWi = prob_tree(vProds) / insideTableKL[0][0][words[i].size() - 1][0] ;

        auto startUpdateTheta = std::chrono::system_clock::now();
        double pTiTiMinis1 = p_ti_ti_minus_1_vec_opt(i);

        updateThetaTime += std::chrono::system_clock::now() - startUpdateTheta;
        vProds.clear();


        sample_parse_tree_kl_vec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size() - 1);
        std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> tp;
        std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> tpaux;
        tpaux = parse_trees_vec[i];
        tp.first = parse_trees_vec[i].first,
        tp.second = vProds;
        parse_trees_vec[i] = tp;
        double pTiLineWi = prob_tree(vProds) / insideTableKL[0][0][words[i].size() - 1][0] ;

        double pTiLineTiMinis1 = p_ti_ti_minus_1_vec_opt(i);
        double funcA = std::min(1.0,(pTiLineTiMinis1*pTiWi)/(pTiTiMinis1*pTiLineWi));
        free_inside_table_kl(insideTableKL, words[i].size());
        p = dist(mt);
        if (funcA >= p) {
            parse_trees_vec[i] = tpaux;
            p_ti_minus_1_plus_frequence(i);
        }
        else{
            countAcceptedTress++;
        }
        remainingTime += std::chrono::system_clock::now() - startRemaining;
        totalTime += std::chrono::system_clock::now() - startIt;
    }
    std::cout << "   iteration: " << iterations  << " Accepted Tress: " << countAcceptedTress << std::endl;
    //print_grammar();
}

void Grammar::Grammar::generate_rules_cnf() {

    for (unsigned long i = 0; i <= context_size.first; i++) {
        for (unsigned long j = 0; j <= context_size.second; j++) {
            std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
            std::vector<std::vector<Symbol::Symbol>> permutationsAll;
            std::vector<std::vector<Symbol::Symbol>> permutationsNonterminals;
            std::vector<Symbol::Symbol> word;
            std::vector<Symbol::Symbol> allSymbols;
            allSymbols.insert(allSymbols.end(), terminals.begin(), terminals.end());
            allSymbols.insert(allSymbols.end(), non_terminals.begin(), non_terminals.end());
            generate_permutation(permutationsTerminals, terminals, i, word, true);
            generate_permutation(permutationsAll, non_terminals, j, word, true);
            generate_permutation(permutationsNonterminals, non_terminals, 2, word, false);
            std::vector<std::vector<Symbol::Symbol>>::iterator itPermutationsTerminals;
            std::vector<std::vector<Symbol::Symbol>>::iterator itPermutationsAll;
            std::vector<std::vector<Symbol::Symbol>>::iterator itPermutationsNonterminals;
            std::vector<Symbol::Symbol>::iterator itNonterminals;
            std::vector<Symbol::Symbol>::iterator itTerminals;
            for( itPermutationsTerminals = permutationsTerminals.begin(); itPermutationsTerminals != permutationsTerminals.end(); itPermutationsTerminals++) {
                for (itNonterminals = non_terminals.begin(); itNonterminals != non_terminals.end(); itNonterminals++) {
                    for (itPermutationsAll = permutationsAll.begin(); itPermutationsAll != permutationsAll.end(); itPermutationsAll++) {
                        std::vector<Symbol::Symbol> leftHandSide;
                        leftHandSide.insert(leftHandSide.end(), (*itPermutationsTerminals).begin(), (*itPermutationsTerminals).end());
                        leftHandSide.push_back(*itNonterminals);
                        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> rulesByLeft;
                        leftHandSide.insert(leftHandSide.end(), (*itPermutationsAll).begin(), (*itPermutationsAll).end());
                        for (itPermutationsNonterminals = permutationsNonterminals.begin(); itPermutationsNonterminals !=  permutationsNonterminals.end(); itPermutationsNonterminals++) {
                            std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rightHandSide;
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsTerminals).begin(), (*itPermutationsTerminals).end());
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsNonterminals).begin(), (*itPermutationsNonterminals).end());
                            rightHandSide.second.first = rightHandSide.second.second = 0.0;
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsAll).begin(), (*itPermutationsAll).end());
                            rulesByLeft.push_back(rightHandSide);
                        }
                        for (itTerminals = terminals.begin(); itTerminals !=  terminals.end(); itTerminals++) {
                            std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightHandSide;
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsTerminals).begin(), (*itPermutationsTerminals).end());
                            Symbol::Symbol aux = (*itTerminals).clone();
                            aux.context = false;
                            rightHandSide.first.push_back(aux);
                            rightHandSide.second.first = rightHandSide.second.second = 0.0;
                            rightHandSide.first.insert(rightHandSide.first.end(), (*itPermutationsAll).begin(), (*itPermutationsAll).end());
                            rulesByLeft.push_back(rightHandSide);
                        }
                        Rule::Rule r = Rule::Rule(leftHandSide, rulesByLeft);
                        r.generate_pior_dirichlet(ALFA);
                        r.update_prob_dirichlet_theta();
                        r.left_context.insert(r.left_context.end(), (*itPermutationsTerminals).begin(), (*itPermutationsTerminals).end());
                        r.right_context.insert(r.right_context.end(), (*itPermutationsAll).begin(), (*itPermutationsAll).end());
                        r.index_1st_non_context = r.left_context.size();
                        rules.push_back(r);
                    }
                }
            }
        }
    }
    std::vector<Rule::Rule>::iterator itRule;
    int j = 0;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        (*itRule).index = j;
        for (unsigned long i =0; i < (*itRule).right.size(); i++)
            (*itRule).rule_frequence.emplace_back(ALFA, 0);
        j++;
    }


}


int Grammar::Grammar::convert_context_to_id(int side, std::vector<Symbol::Symbol> context) {
    size_t nSymbols;
    int id = 0;
    if (context.empty())
        return id;
    if (side == 0)
        nSymbols = n_terminals;
    else
        nSymbols = n_non_terminals;
    std::vector<Symbol::Symbol>::iterator itSymbol;
    for (unsigned long i = 0; i < context.size(); i++)
        id += (int) pow(nSymbols, i);
    int position = 0;
    for (itSymbol = context.begin(); itSymbol != context.end(); itSymbol++) {
        id += static_cast<int>((*itSymbol).id * pow(nSymbols, position));
        position++;
    }
    return id;
}

void Grammar::Grammar::sample_parse_tree_kl(
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> &vr,
        Rule::Rule r, const std::string& w, double ****inside_table, size_t i, size_t k) {
    size_t jRange = k - i;
    size_t rightRange = n_non_terminals * n_non_terminals;

    Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
    std::vector<Symbol::Symbol>::iterator itLeft;
    std::vector<Symbol::Symbol> leftContext;
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
    std::vector<Symbol::Symbol> rightContext;
    while (itLeft != r.left.end() ) {
        rightContext.push_back((*itLeft));
        itLeft++;
    }
    //double sumJBC = 0.0;
    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    std::vector<Symbol::Symbol> rightContextLine;
                    itSymbol++;
                    if (context_size.second > 0) {
                        rightContextLine.push_back((*itSymbol));
                        rightContextLine.front().context = true;
                        if (!rightContext.empty()) {
                            rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                            rightContext.pop_back();
                        }
                    }
                    itSymbol--;
                    for (unsigned long l = 0; l < jRange; l++) {
                        size_t indexB = (*itSymbol).id * n_non_terminals * jRange;
                        double pBij = inside_table[context_amount.first * (*itSymbol).id +
                                convert_context_to_id(0, leftContext)]
                                [i][i + l][convert_context_to_id(1, rightContextLine)];
                        itSymbol++;
                        double pCjk = inside_table[context_amount.first * (*itSymbol).id +
                                convert_context_to_id(0, leftContext)]
                                [i + l + 1][k][convert_context_to_id(1, rightContext)];
                        size_t indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) / inside_table
                                [context_amount.first * nonterminal.id + convert_context_to_id(0, leftContext)]
                                [i][k][convert_context_to_id(1, rightContext)];
                        itSymbol--;
                        /*std::cout << "SumJBC = " <<sumJBC << " ";

                        std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first <<
                        " * " << pBij << " ("<< "IT["<<
                                      left_context.size()*(*itSymbol).id+ convertContextToID(0,left_context)<<"]["<<i<<"]["<<i+l<<"]["<<
                                      convertContextToID(1,right_context)<<"])"
                        << " * " << pCjk
                                << " ("<< "IT["<<
                                left_context.size()*(*itSymbol).id+ convertContextToID(0,left_context)<<"]["<<i+l+1<<"]["<<k<<"]["<<
                                convertContextToID(1,right_context)<<"])"
                        << " / " << inside_table
                        [left_context.size()*nonterminal.id+ convertContextToID(0,left_context)]
                        [i][k][convertContextToID(1,right_context)] << " IT["<<
                        left_context.size()*nonterminal.id+ convert_context_to_id(0,left_context)<<"]["<<i<<"]["<<k<<"]["<<
                        convertContextToID(1,right_context)<<"]"<<std::endl;*/


                    }
                    break;
                }
            }
        }


        unsigned long l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //std::cout << "Sensitive: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        size_t j = 0;

        delete[]jbc;
        Symbol::Symbol rightB = Symbol::Symbol("", 0,false);
        Symbol::Symbol rightC = Symbol::Symbol("", 0,false);

        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    if ((*itSymbol).id == l / (n_non_terminals * jRange)) {
                        itSymbol++;
                        if(itSymbol == (*itRight).first.end())
                            break;
                        if((*itSymbol).id == (l%(n_non_terminals * jRange)) / jRange) {
                            itSymbol--;
                            rightProduction.first.push_back((*itSymbol));
                            rightB = (*itSymbol).clone();
                            j = l - n_non_terminals * jRange * (*itSymbol).id;
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
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        apply_production(production);
        std::vector<Symbol::Symbol> leftContext2;
        std::vector<Symbol::Symbol> rightContext2;
        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();

        size_t lc , rc;
        std::uniform_int_distribution<size_t> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<size_t> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = ((int) rand()%(leftContext2.size()+1));
        //rc = (int) rand()%(rightContext2.size()+1);
        for (size_t i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (size_t i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::vector<Symbol::Symbol> lhs;
        lhs.insert(lhs.end(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightB);
        lhs.insert(lhs.end(), rightContext2.begin(), rightContext2.end());
        Rule::Rule nextRule = find_rule_by_lhs(lhs);
        sample_parse_tree_kl(vr, nextRule, w, inside_table, i, i + j);

        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();


        lc = distIntL(mt);
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (unsigned int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (unsigned int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();

        lhs.clear();
        lhs.insert(lhs.begin(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightC);
        lhs.insert(lhs.begin(), rightContext2.begin(), rightContext2.end());
        sample_parse_tree_kl(vr, find_rule_by_lhs(lhs), w, inside_table, i + j + 1, k);

    } else {
        std::vector<Symbol::Symbol> leftContext2;
        std::vector<Symbol::Symbol> rightContext2;
        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();
        size_t lc , rc;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<size_t> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<size_t> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (size_t i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (size_t i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();

        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w.substr(i,jRange+1);
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
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
        //rightProduction.first.insert(rightProduction.first.end(), right_context.begin(), right_context.end());
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        apply_production(production);
        return;
    }
}

void Grammar::Grammar::apply_production(
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> prod_a) {
    std::vector<Symbol::Symbol>::iterator itLeft;
    std::vector<Symbol::Symbol> newProd;
    bool appliedProd = false;
    for (itLeft = prod_a.first.begin(); itLeft != prod_a.first.end(); itLeft++) {
        if (!itLeft->terminal) {
            std::vector<Symbol::Symbol>::iterator itActualProd;

            for (itActualProd = actual_production.begin(); itActualProd != actual_production.end(); itActualProd++) {
                if ((*itActualProd).equal_symbol(*itLeft) && !appliedProd) {
                    std::vector<Symbol::Symbol>::iterator  itRight;
                    for (itRight = prod_a.second.first.begin(); itRight != prod_a.second.first.end(); itRight++) {
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
    actual_production = newProd;
}

Rule::Rule Grammar::Grammar::find_rule_by_lhs(std::vector<Symbol::Symbol> lhs) {
    std::vector<Rule::Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<Symbol::Symbol>::iterator itLHS;
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

void Grammar::Grammar::get_actual_context(std::vector<Symbol::Symbol> &leftContext, std::vector<Symbol::Symbol> &rightContext) {
    leftContext.clear();
    rightContext.clear();
    std::vector<Symbol::Symbol>::iterator itLeft;
    Symbol::Symbol aux = Symbol::Symbol("",0,true);
    for (itLeft = actual_production.begin(); itLeft != actual_production.end(); itLeft++) {
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
    while (itLeft != actual_production.end() ) {
        aux = *itLeft;
        aux.context = true;
        aux.terminal = false;
        rightContext.push_back(aux);
        itLeft++;
    }
}

void Grammar::Grammar::sample_parse_tree_vec(
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> & vr,
        Rule::Rule r, std::vector<Symbol::Symbol> w, double ***p_double, size_t i, size_t k) {
    size_t jRange = k - i;
    size_t rightRange = n_non_terminals * n_non_terminals;

    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    for (unsigned int l = 0; l < jRange; l++) {
                        size_t indexB = (*itSymbol).id * n_non_terminals * jRange;
                        double pBij = p_double[(*itSymbol).id][i][i + l];
                        itSymbol++;
                        double pCjk = p_double[(*itSymbol).id][i + l + 1][k];
                        size_t indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) / p_double[(r.left[0].id)][i][k];
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
        unsigned int l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //std::cout << "Free: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        size_t j = 0;
        delete[]jbc;
        //std::vector<Symbol::Symbol> rightProduction;
        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    if ((*itSymbol).id == l / (n_non_terminals * jRange)) {
                        itSymbol++;
                        if(itSymbol == (*itRight).first.end())
                            break;
                        if((*itSymbol).id == (l%(n_non_terminals * jRange)) / jRange) {
                            itSymbol--;
                            rightProduction.first.push_back((*itSymbol));
                            j = l - n_non_terminals * jRange * (*itSymbol).id;
                            itSymbol++;
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
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        sample_parse_tree_vec(vr, rules[production.second.first[0].id], w, p_double, i, i + j);
        sample_parse_tree_vec(vr, rules[production.second.first[1].id], w, p_double, i + j + 1, k);

    } else {
        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w[i].name;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            //TO DO:: Dar um jeito de não hardecodá o left[0]
            if ((*itRight).first[0].name == s) {
                rightProduction.first.push_back((*itRight).first[0]);
                rightProduction.second.first = (*itRight).second.first;
                rightProduction.second.second = (*itRight).second.second;
            }
        }
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        return;
    }

}

double ***Grammar::Grammar::cyk_prob_vec(std::vector<Symbol::Symbol> w) {
    //std::cout <<"IT for " << w << std::endl;
    auto ***p = new double**[non_terminals.size()];
    for (unsigned long i = 0; i < non_terminals.size(); i++) {
        p[i] = new double *[w.size()];
        for (unsigned long j = 0; j < w.size(); j++)
            p[i][j] = new double[w.size()];
    }

    for (unsigned long i = 0; i < non_terminals.size(); i++)
        for (unsigned long j = 0; j < w.size(); j++)
            for (unsigned long k = 0; k < w.size(); k++)
                p[i][j][k] = 0.0;

    std::vector<Rule::Rule>::iterator itRule;
    for (unsigned long i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol::Symbol>::iterator itLeft;
            Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context)
                    nonterminal = (*itLeft).clone();
                break;
            }

            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol::Symbol>::iterator itRightS;
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
    for (unsigned long i = 1; i < w.size(); i++) {
        for (unsigned long j = 0; j < w.size()-i; j++) {
            for (unsigned long k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol::Symbol>::iterator itLeft;
                    Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
                    for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                        if (!itLeft->terminal && !itLeft->context)
                            nonterminal = (*itLeft).clone();
                        break;
                    }
                    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol::Symbol>::iterator itRightS;
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

    //print_inside_table(p,w.size());
    return p;
}

vector<vector<vector<double>>> Grammar::Grammar::outside_table(const std::vector<Symbol::Symbol>& w, double *** inside_table) {
    vector<vector<vector<double>>>  outsideTable;
    for (const auto& nt: non_terminals) {
        outsideTable.emplace_back();
        for (unsigned long i = 0; i < w.size(); i++) {
            outsideTable[nt.id].push_back(vector<double>());
            for (unsigned long j = 0; j < w.size(); j++) {
                outsideTable[nt.id][i].push_back(0.0);
            }
        }
    }
    outsideTable[0][0][w.size()-1] = 1;
    //print_outside_table(outside_table);

    for (int i = (int) w.size()-1; i >=0; i--) {
        for (int j = 0; j+i< (int) w.size(); j++) {
            for (const auto& nt: non_terminals) {
                //cout << "B_"<<nt.id<<"_"<<j<<"_"<<j+i<<endl;

                for (auto r: rules) {
                    for (auto right: r.right) {
                        if (right.first[0].terminal)
                            break;
                        for (int k=0; k < j; k++) {

                            if (right.first[1].id ==nt.id) {
                                outsideTable[nt.id][j][j+i] += outsideTable[r.left[0].id][k][j+i] * inside_table[ right.first[0].id][k][j - 1] * right.second.first;
                                /*cout << "B_" << r.left[0].id << "_" << k << "_" << j + i << " * " << "I_" << right.first[0].id << "_" << k << "_" << j - 1 << " = "
                                     <<outsideTable[r.left[0].id][k][j+i] << " * " << inside_table[ right.first[0].id][k][j-1] << " * " << right.second.first
                                     << " = " << outside_table[nt.id][j][j+i] << endl;*/
                            }                    }
                        for (int k= (int) w.size()-1; k > j+i; k--) {
                            if (right.first[0].id ==nt.id) {
                                outsideTable[nt.id][j][j+i] += outsideTable[r.left[0].id][j][k] * inside_table[ right.first[1].id][j + i + 1][k] * right.second.first;
                                /*cout << "B_" <<r.left[0].id<<"_"<<j<<"_"<<k<< " * " << "I_"<<right.first[1].id<<"_"<<j+i+1<<"_"<<k << " = "
                                     << outsideTable[r.left[0].id][j][k] << " * " <<inside_table[ right.first[1].id][j+i+1][k] << " * " << right.second.first
                                     << " = " << outside_table[nt.id][j][j+i] <<endl;*/
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
                for (int l = 0; l < n_non_terminals*n_non_terminals; l++) {
                    for (int s = j; s <=k; s++) {
                        outsideTable[r.right[l].first[0].id][j][s] += r.right[l].second.first *
                                outsideTable[r.left[0].id][j][k] * inside_table[r.right[l].first[1].id][s][k-1];

                        cout << "O["<<r.right[l].first[0].id<<"]["<<j<<"]["<<s<<"] += O[" << r.left[0].id<<"]["<<j<<"]["<<k<<"] * I["
                        <<r.right[l].first[1].id<<"]["<<s<<"]["<<k-1<<"]\t\t" ;



                        outsideTable[r.right[l].first[1].id][s][k] += r.right[l].second.first *
                                outside_table[r.left[0].id][j][k] * inside_table[r.right[l].first[0].id][j+1][s];

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
                for (int l = 0; l < n_non_terminals*n_non_terminals; l++) {
                    for (int k = q+1; k <w.size();k++) {
                        outsideTable[r.right[l].first[0].id][p][q] += outside_table[r.left[0].id][p][k] *
                                inside_table[r.right[l].first[1].id][q+1][k]*r.right[l].second.first;
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
                    for (int l = 0; l < n_non_terminals*n_non_terminals; l++) {
                        cout << "i:"<<i<<" j:"<<j<<" k:"<<k<<" - ";
                        outside_table[r.right[l].first[0].id][j][j+k] +=
                                outsideTable[r.left[0].id][j][j+i] * inside_table[r.right[l].first[1].id][j+k][j+i] * r.right[l].second.first;

                        cout << "OT["<<r.right[l].first[0].id<<"]["<<j<<"]["<<j+k<<"] += OT["<<r.left[0].id<<"]["<<j<<"]["<<j+i<<"] * IT["
                             << r.right[l].first[1].id<<"]["<<j+k<<"]["<<j+i<<"] * p"<<r.left[0].id<<"->"<<r.right[l].first[0].id<<r.right[l].first[1].id << " = " << outsideTable[r.left[0].id][j][j+i] << "*"
                             << inside_table[r.right[l].first[1].id][j+k][j+i] <<"*" << r.right[l].second.first << " = " << outsideTable[r.right[l].first[0].id][j][j+k];

                        outsideTable[r.right[l].first[1].id][j+k][j+i] +=
                                outsideTable[r.left[0].id][j][j+i] * inside_table[r.right[l].first[0].id][j][j+k] * r.right[l].second.first;

                        cout << "\tOT["<<r.right[l].first[1].id<<"]["<<j+k<<"]["<<j+i<<"] += OT["<<r.left[0].id<<"]["<<j<<"]["<<j+i<<"] * IT["
                             <<r.right[l].first[0].id<<"]["<<j<<"]["<<j+k<<"] * p"<<r.left[0].id<<"->"<<r.right[l].first[0].id<<r.right[l].first[1].id << " = " << outsideTable[r.left[0].id][j][j+i] << "*"
                             << inside_table[r.right[l].first[0].id][j][j+k] <<"*" << r.right[l].second.first << " = " << outsideTable[r.right[l].first[0].id][j+k][j+i] << endl;
                    }
                }
            }
        }
    }*/
    return outsideTable;
}

void Grammar::Grammar::calculate_new_theta_vec(const std::vector<Symbol::Symbol>& w) {
    std::vector<Rule::Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        std::vector<int> nonterminalCounts;
        double nonterminalTotal = 0;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production;
            production.first = (*itRule).left;
            production.second.first = (*itRight).first;
            int productionCount = calculate_producton_counts_vec(production, w);
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
        //(*itRule).print_rule();

    }
}

int Grammar::Grammar::calculate_producton_counts_vec(
        const std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>& pair,
        const std::vector<Symbol::Symbol>& w) {
    int productionCount = 0;

    std::vector<std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>>>::iterator itParseTreeVec;
    for (itParseTreeVec = parse_trees_vec.begin(); itParseTreeVec != parse_trees_vec.end(); itParseTreeVec++) {
        if (!equal_word(w, (*itParseTreeVec).first)) {
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>::iterator itProduction;
            for (itProduction = (*itParseTreeVec).second.begin(); itProduction != (*itParseTreeVec).second.end(); itProduction++) {
                if(equal_productions((*itProduction), pair)) {
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

bool Grammar::Grammar::equal_word(std::vector<Symbol::Symbol> w1, std::vector<Symbol::Symbol> w2) {
    std::vector<Symbol::Symbol>::iterator itW1;
    auto itW2 = w2.begin();
    for (itW1 = w1.begin(); itW1 != w1.end(); itW1++) {
        if(itW2 == w2.end())
            return false;
        if (!(*itW1).equal_symbol((*itW2)))
            return false;
        itW2++;
    }
    if(itW2!= w2.end())
        return false;
    else
        return true;


}

double Grammar::Grammar::p_ti_ti_minus_1_vec(const std::vector<Symbol::Symbol>& w) {
    std::vector<Rule::Rule>::iterator itRule;
    double result = 1.0;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<Symbol::Symbol> a;
        result *= c_constant(calculate_rule_frequence_vec((*itRule), a)) /
                  c_constant(calculate_rule_frequence_vec((*itRule), w));
        //std::cout << "   PTia = " << c_constant(calculateRuleFrequenceVec((*itRule),  a)) << " Pt-ia = " << c_constant(calculate_rule_frequence_vec((*itRule), w)) << std::endl;
        //Talvez essa passagem de parametro esteja errada

    }
    return result;
}

double Grammar::Grammar::p_ti_ti_minus_1_vec_opt(int i) {
    std::vector<Rule::Rule>::iterator itRule;
    std::vector<std::vector<std::pair<double, int>>> vecPairs;
    vecPairs.reserve(rules.size());
    for (itRule = rules.begin();itRule < rules.end(); itRule++) {
        vecPairs.push_back(itRule->rule_frequence);
    }
    p_ti_minus_1_frequence(i);
    double result = 1.0;
    auto vecPairsIt = vecPairs.begin();
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<Symbol::Symbol> a;
        result *= c_constant(*vecPairsIt) / c_constant(itRule->rule_frequence);
        vecPairsIt++;
    }
    return result;
}

void Grammar::Grammar::p_ti_minus_1_frequence(int i) {
    if (i == -1)
        return;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>::iterator pTreeIt;
    for (pTreeIt = parse_trees_vec[i].second.begin(); pTreeIt != parse_trees_vec[i].second.end(); pTreeIt++) {
        Rule::Rule r = find_rule_by_lhs(pTreeIt->first);
        if (!(*pTreeIt).second.first[r.index_1st_non_context].terminal)
            rules[r.index].rule_frequence[(*pTreeIt).second.first[r.index_1st_non_context].id * n_non_terminals +
                                          (*pTreeIt).second.first[r.index_1st_non_context + 1].id].second--;
        else
            rules[r.index].rule_frequence[n_non_terminals * n_non_terminals + (*pTreeIt).second.first[r.index_1st_non_context].id].second--;
    }
}

void Grammar::Grammar::p_ti_minus_1_plus_frequence(int i) {
    if (i == -1)
        return;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>::iterator pTreeIt;
    for (pTreeIt = parse_trees_vec[i].second.begin(); pTreeIt != parse_trees_vec[i].second.end(); pTreeIt++) {
        Rule::Rule r = find_rule_by_lhs(pTreeIt->first);
        if (!(*pTreeIt).second.first[r.index_1st_non_context].terminal)
            rules[r.index].rule_frequence[(*pTreeIt).second.first[r.index_1st_non_context].id * n_non_terminals +
                                          (*pTreeIt).second.first[r.index_1st_non_context + 1].id].second++;
        else
            rules[r.index].rule_frequence[n_non_terminals * n_non_terminals + (*pTreeIt).second.first[r.index_1st_non_context].id].second++;
    }
}

std::vector<std::pair<double, int>> Grammar::Grammar::calculate_rule_frequence_vec(Rule::Rule &rule, const std::vector<Symbol::Symbol>& w) {
    std::vector<std::pair<double, int>> ruleFrequence;
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = rule.right.begin(); itRight != rule.right.end(); itRight++) {
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production;
        production.first = rule.left;
        production.second.first = (*itRight).first;
        std::pair<double, int> productionFrequence;
        productionFrequence.first = (*itRight).second.second;
        productionFrequence.second = calculate_producton_counts_vec(production, w);
        ruleFrequence.push_back(productionFrequence);
    }
    return ruleFrequence;
}



double ****Grammar::Grammar::cyk_prob_kl_vec(std::vector<Symbol::Symbol> w) {
    size_t sizeLeftContext = 0;
    size_t sizeRightContext = 0;
    for (unsigned long i = 0; i <= context_size.first; i++) {
        std::vector<Symbol::Symbol> word;
        std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
        generate_permutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();

    }
    context_amount.first = sizeLeftContext;
    for (size_t i = 0; i <= context_size.second; i++) {
        std::vector<Symbol::Symbol> word;
        std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
        generate_permutation(permutationsTerminals, non_terminals, i, word, true);
//        rightContexts.insert(rightContexts.end(), permutationsTerminals.begin(), permutationsTerminals.end());
        sizeRightContext += permutationsTerminals.size();
    }
    context_amount.second = sizeRightContext;
    auto ****p = new double***[non_terminals.size() * sizeLeftContext];
    for (unsigned long i = 0; i < non_terminals.size() * sizeLeftContext; i++) {
        p[i] = new double **[w.size()];
        for (unsigned long j = 0; j < w.size(); j++) {
            p[i][j] = new double *[w.size()];
            for (unsigned long k = 0; k < w.size(); k++)
                p[i][j][k] = new double[sizeRightContext];
        }
    }

    for (unsigned long i = 0; i < non_terminals.size() * sizeLeftContext; i++)
        for (unsigned long j = 0; j < w.size(); j++)
            for (unsigned long k = 0; k < w.size(); k++)
                for (size_t l = 0; l < sizeRightContext; l++)
                    p[i][j][k][l] = 0.0;

    std::vector<Rule::Rule>::iterator itRule;
    for (unsigned long i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol::Symbol>::iterator itLeft;
            Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
            std::vector<Symbol::Symbol> leftContext;
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
            std::vector<Symbol::Symbol> rightContext;
            while (itLeft != (*itRule).left.end() ) {
                rightContext.push_back((*itLeft));
                itLeft++;
            }

            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol::Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                        if(itRightS->name == w[i].name) {
                            p[sizeLeftContext*nonterminal.id + convert_context_to_id(0, leftContext)][i][i][convert_context_to_id(
                                    1, rightContext)] = itRight->second.first;
                            //std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,left_context)<<"]["<<i<<"]["<<i<<"]["<< convert_context_to_id(1,right_context)<<"] = ProbR"
                            /* <<(*itRight).second.first << std::endl;*/
                            break;
                        }
                    }
                }
            }

        }
    }
    for (unsigned long i = 1; i < w.size(); i++) {
        for (unsigned long j = 0; j < w.size()-i; j++) {
            for (unsigned long k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol::Symbol>::iterator itLeft;
                    Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
                    std::vector<Symbol::Symbol> leftContext;
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
                    std::vector<Symbol::Symbol> rightContext;
                    while (itLeft != (*itRule).left.end() ) {
                        rightContext.push_back((*itLeft));
                        itLeft++;
                    }

                    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol::Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (!itRightS->terminal && !itRightS->context) {
                                itRightS++;
                                std::vector<Symbol::Symbol> rightContextLine;
                                if (context_size.second > 0) {
                                    rightContextLine.push_back((*itRightS));
                                    rightContextLine.front().context = true;
                                    if (!rightContext.empty()) {
                                        rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                                        rightContextLine.pop_back();
                                    }
                                }
                                int leftContextId = convert_context_to_id(0, leftContext);
                                int rightContextID = convert_context_to_id(1, rightContext);
                                itRightS--;
                                double bInside = p[sizeLeftContext*(*itRightS).id + leftContextId][j][j+k][convert_context_to_id(
                                        1, rightContextLine)];
                                itRightS++;
                                double cInside = p[sizeLeftContext*(*itRightS).id + leftContextId][j+k+1][j+i][rightContextID];
                                p[sizeLeftContext*nonterminal.id+ leftContextId][j][i+j][rightContextID] =
                                        p[sizeLeftContext*nonterminal.id+ leftContextId][j][i+j][rightContextID] + bInside*cInside*(*itRight).second.first;

                                /*std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,left_context)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,right_context)<<"] = "
                                         << "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,left_context)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,right_context)<<
                                         "] + IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,left_context)<<"]["<<j<<"]["<<j+k<<"]["<< convertContextToID(1,right_context)<<
                                         "] * IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,left_context)<<"]["<<j+k+1<<"]["<<j+i<<"]["<< convertContextToID(1,right_context) << "] * " <<  "ProbR = "
                                         <<p[sizeLeftContext*nonterminal.id+ convertContextToID(0,left_context)][j][i+j][convert_context_to_id(1,right_context)]<< std::endl;*/
                                //std::cout << p[sizeLeftContext*nonterminal.id+ convertContextToID(0,left_context)][j][i+j][convert_context_to_id(1,right_context)]  << " + " << bInside << " * " << cInside << " * " << (*itRight).second.first << std::endl;
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

void Grammar::Grammar::sample_parse_tree_kl_vec(
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> &vr,
        Rule::Rule r, std::vector<Symbol::Symbol> w, double ****inside_table, size_t i, size_t k) {
    size_t jRange = k - i;
    size_t rightRange = n_non_terminals * n_non_terminals;

    Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
    std::vector<Symbol::Symbol>::iterator itLeft;
    std::vector<Symbol::Symbol> leftContext;
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
    std::vector<Symbol::Symbol> rightContext;
    while (itLeft != r.left.end() ) {
        rightContext.push_back((*itLeft));
        itLeft++;
    }
    //double sumJBC = 0.0;
    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    std::vector<Symbol::Symbol> rightContextLine;
                    itSymbol++;
                    if (context_size.second > 0) {
                        rightContextLine.push_back((*itSymbol));
                        rightContextLine.front().context = true;
                        if (!rightContext.empty()) {
                            rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                            rightContextLine.pop_back();
                        }
                    }
                    itSymbol--;
                    for (unsigned int l = 0; l < jRange; l++) {
                        size_t indexB = (*itSymbol).id * n_non_terminals * jRange;
                        double pBij = inside_table[context_amount.first * (*itSymbol).id +
                                convert_context_to_id(0, leftContext)]
                        [i][i + l][convert_context_to_id(1, rightContextLine)];
                        itSymbol++;
                        double pCjk = inside_table[context_amount.first * (*itSymbol).id +
                                convert_context_to_id(0, leftContext)]
                        [i + l + 1][k][convert_context_to_id(1, rightContext)];
                        size_t indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) / inside_table
                        [context_amount.first * nonterminal.id + convert_context_to_id(0, leftContext)]
                        [i][k][convert_context_to_id(1, rightContext)];
                        itSymbol--;
                        /*std::cout << "SumJBC = " <<sumJBC << " ";

                        std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first <<
                        " * " << pBij << " ("<< "IT["<<
                                      left_context.size()*(*itSymbol).id+ convertContextToID(0,left_context)<<"]["<<i<<"]["<<i+l<<"]["<<
                                      convertContextToID(1,right_context)<<"])"
                        << " * " << pCjk
                                << " ("<< "IT["<<
                                left_context.size()*(*itSymbol).id+ convertContextToID(0,left_context)<<"]["<<i+l+1<<"]["<<k<<"]["<<
                                convertContextToID(1,right_context)<<"])"
                        << " / " << inside_table
                        [left_context.size()*nonterminal.id+ convertContextToID(0,left_context)]
                        [i][k][convertContextToID(1,right_context)] << " IT["<<
                        left_context.size()*nonterminal.id+ convertContextToID(0,left_context)<<"]["<<i<<"]["<<k<<"]["<<
                        convert_context_to_id(1,right_context)<<"]"<<std::endl;*/


                    }
                    break;
                }
            }
        }


        unsigned int l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //std::cout << "Sensitive: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        size_t j = 0;

        delete[]jbc;
        Symbol::Symbol rightB = Symbol::Symbol("", 0,false);
        Symbol::Symbol rightC = Symbol::Symbol("", 0,false);

        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    if ((*itSymbol).id == l / (n_non_terminals * jRange)) {
                        itSymbol++;
                        if(itSymbol == (*itRight).first.end() || (*itSymbol).context)
                            break;
                        if((*itSymbol).id == (l%(n_non_terminals * jRange)) / jRange) {
                            itSymbol--;
                            rightProduction.first.push_back((*itSymbol));
                            rightB = (*itSymbol).clone();
                            j = l - n_non_terminals * jRange * (*itSymbol).id;
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

        rules[r.index].rule_frequence[rightB.id * n_non_terminals + rightC.id].second++;
        rightProduction.first.insert(rightProduction.first.end(), rightContext.begin(), rightContext.end());
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        apply_production(production);
        std::vector<Symbol::Symbol> leftContext2;
        std::vector<Symbol::Symbol> rightContext2;
        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();

        size_t lc , rc;
        std::uniform_int_distribution<size_t> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<size_t> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = ((int) rand()%(leftContext2.size()+1));
        //rc = (int) rand()%(rightContext2.size()+1);
        for (size_t i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (size_t i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::vector<Symbol::Symbol> lhs;
        lhs.insert(lhs.end(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightB);
        lhs.insert(lhs.end(), rightContext2.begin(), rightContext2.end());
        Rule::Rule nextRule = find_rule_by_lhs(lhs);
        sample_parse_tree_kl_vec(vr, nextRule, w, inside_table, i, i + j);

        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();


        lc = distIntL(mt);
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (unsigned int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (unsigned int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();

        lhs.clear();
        lhs.insert(lhs.begin(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightC);
        lhs.insert(lhs.end(), rightContext2.begin(), rightContext2.end());
        sample_parse_tree_kl_vec(vr, find_rule_by_lhs(lhs), w, inside_table, i + j + 1, k);

    } else {
        std::vector<Symbol::Symbol> leftContext2;
        std::vector<Symbol::Symbol> rightContext2;
        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();
        size_t lc , rc;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<size_t> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<size_t> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (size_t i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (size_t i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w[i].name;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            bool flagContext = false;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if ((*itSymbol).name == s && !(*itSymbol).context ) {
                    rightProduction.first.push_back((*itSymbol));
                    rightProduction.second.first = (*itRight).second.first;
                    rightProduction.second.second = (*itRight).second.second;
                    flagContext = true;
                    rules[r.index].rule_frequence[n_non_terminals * n_non_terminals + (*itSymbol).id].second++;
                } else if (flagContext)
                    rightProduction.first.push_back((*itSymbol));
            }

        }
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        apply_production(production);
        return;
    }
}

std::pair<double,double> Grammar::Grammar::perplexity(const std::vector<std::vector<Symbol::Symbol>>& test_data) {
    std::pair<double,double> pXTheta = std::make_pair(1.0, 0.0);
    for (const auto& w: test_data) {
        double ***iTable  = cyk_prob_vec(w);
        double pX = iTable[0][0][w.size()-1] ;
        double pXN = pX;
        double ***pN = cyk_prob_n_vec(w);
        pXN /= pN[0][0][w.size()-1];
        free_inside_table(pN, w.size());
        free_inside_table(iTable, w.size());
        pXTheta.first += log(pX);
        pXTheta.second += log(pXN);

    }
    return std::make_pair(exp(-(pXTheta.first / test_data.size())), exp(-(pXTheta.second / test_data.size())));
}

std::pair<double,double> Grammar::Grammar::perplexity_kl(const std::vector<std::vector<Symbol::Symbol>>& test_data) {
    std::pair<double,double> pXTheta = std::make_pair(1.0, 0.0);
    for (const auto& w: test_data) {
        double ****iTable  = cyk_prob_kl_vec(w);
        double pX = iTable[0][0][w.size()-1][0];
        double pXN = pX;
        double ****pN = cyk_prob_kln_vec(w);
        pXN /= pN[0][0][w.size()-1][0];
        free_inside_table_kl(pN, w.size());
        free_inside_table_kl(iTable, w.size());
        pXTheta.first += log(pX);
        pXTheta.second += log(pXN);

    }
    return std::make_pair(exp(-(pXTheta.first / test_data.size())), exp(-(pXTheta.second / test_data.size())));
}

double ***Grammar::Grammar::cyk_prob_n_vec(const std::vector<Symbol::Symbol>& w) {
    auto ***p = new double**[non_terminals.size()];
    for (unsigned long i = 0; i < non_terminals.size(); i++) {
        p[i] = new double *[w.size()];
        for (unsigned long j = 0; j < w.size(); j++)
            p[i][j] = new double[w.size()];
    }

    for (unsigned long i = 0; i < non_terminals.size(); i++)
        for (unsigned long j = 0; j < w.size(); j++)
            for (unsigned long k = 0; k < w.size(); k++)
                p[i][j][k] = 0.0;

    std::vector<Rule::Rule>::iterator itRule;
    for (unsigned long i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol::Symbol>::iterator itLeft;
            Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context)
                    nonterminal = (*itLeft).clone();
                break;
            }

            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol::Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                        p[nonterminal.id][i][i] += itRight->second.first;
                        //std::cout<< "IT["<<nonterminal.id<<"]["<<i<<"]["<<i<<"] = " << itRight->second.first << " probR" << std::endl;
                    }
                }
            }

        }
    }
    for (unsigned long i = 1; i < w.size(); i++) {
        for (unsigned long j = 0; j < w.size()-i; j++) {
            for (unsigned long k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol::Symbol>::iterator itLeft;
                    Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
                    for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                        if (!itLeft->terminal && !itLeft->context)
                            nonterminal = (*itLeft).clone();
                        break;
                    }
                    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol::Symbol>::iterator itRightS;
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

    //print_inside_table(p,w.size());
    return p;
}

double ****Grammar::Grammar::cyk_prob_kln_vec(const std::vector<Symbol::Symbol>& w) {
    size_t sizeLeftContext = 0;
    size_t sizeRightContext = 0;
    for (unsigned long i = 0; i <= context_size.first; i++) {
        std::vector<Symbol::Symbol> word;
        std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
        generate_permutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();

    }
    context_amount.first = sizeLeftContext;
    for (unsigned long i = 0; i <= context_size.second; i++) {
        std::vector<Symbol::Symbol> word;
        std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
        generate_permutation(permutationsTerminals, non_terminals, i, word, true);
//        rightContexts.insert(rightContexts.end(), permutationsTerminals.begin(), permutationsTerminals.end());
        sizeRightContext += permutationsTerminals.size();
    }
    context_amount.second = sizeRightContext;
    auto ****p = new double***[non_terminals.size() * sizeLeftContext];
    for (unsigned long i = 0; i < non_terminals.size() * sizeLeftContext; i++) {
        p[i] = new double **[w.size()];
        for (unsigned long j = 0; j < w.size(); j++) {
            p[i][j] = new double *[w.size()];
            for (unsigned long k = 0; k < w.size(); k++)
                p[i][j][k] = new double[sizeRightContext];
        }
    }

    for (unsigned long i = 0; i < non_terminals.size() * sizeLeftContext; i++)
        for (unsigned long j = 0; j < w.size(); j++)
            for (unsigned long k = 0; k < w.size(); k++)
                for (size_t l = 0; l < sizeRightContext; l++)
                    p[i][j][k][l] = 0.0;

    std::vector<Rule::Rule>::iterator itRule;
    for (unsigned long i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol::Symbol>::iterator itLeft;
            Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
            std::vector<Symbol::Symbol> leftContext;
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
            std::vector<Symbol::Symbol> rightContext;
            while (itLeft != (*itRule).left.end() ) {
                rightContext.push_back((*itLeft));
                itLeft++;
            }

            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
            for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                std::vector<Symbol::Symbol>::iterator itRightS;
                for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                    if (itRightS->terminal && !itRightS->context) {
                            p[sizeLeftContext*nonterminal.id + convert_context_to_id(0, leftContext)][i][i][convert_context_to_id(
                                    1, rightContext)] += itRight->second.first;
                            //std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,left_context)<<"]["<<i<<"]["<<i<<"]["<< convert_context_to_id(1,right_context)<<"] = ProbR"
                            /* <<(*itRight).second.first << std::endl;*/
                    }
                }
            }

        }
    }
    for (unsigned long i = 1; i < w.size(); i++) {
        for (unsigned long j = 0; j < w.size()-i; j++) {
            for (unsigned long k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {

                    std::vector<Symbol::Symbol>::iterator itLeft;
                    Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
                    std::vector<Symbol::Symbol> leftContext;
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
                    std::vector<Symbol::Symbol> rightContext;
                    while (itLeft != (*itRule).left.end() ) {
                        rightContext.push_back((*itLeft));
                        itLeft++;
                    }

                    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                        std::vector<Symbol::Symbol>::iterator itRightS;
                        for (itRightS = (*itRight).first.begin(); itRightS != (*itRight).first.end(); itRightS++) {
                            if (!itRightS->terminal && !itRightS->context) {
                                itRightS++;
                                std::vector<Symbol::Symbol> rightContextLine;
                                if (context_size.second > 0) {
                                    rightContextLine.push_back((*itRightS));
                                    rightContextLine.front().context = true;
                                    if (!rightContext.empty()) {
                                        rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                                        rightContext.pop_back();
                                    }
                                }
                                itRightS--;
                                int leftContextID = convert_context_to_id(0, leftContext);
                                int rightContextID = convert_context_to_id(1, rightContext);
                                double bInside = p[sizeLeftContext*(*itRightS).id + leftContextID][j][j+k][convert_context_to_id(
                                        1, rightContextLine)];
                                itRightS++;
                                double cInside = p[sizeLeftContext*(*itRightS).id + leftContextID][j+k+1][j+i][rightContextID];
                                p[sizeLeftContext*nonterminal.id+ leftContextID][j][i+j][rightContextID] =
                                        p[sizeLeftContext*nonterminal.id+ leftContextID][j][i+j][rightContextID] + bInside*cInside*(*itRight).second.first;

                                /*std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convert_context_to_id(0,left_context)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,right_context)<<"] = "
                                         << "IT["<<sizeLeftContext*nonterminal.id + convertContextToID(0,left_context)<<"]["<<j<<"]["<<i+j<<"]["<< convertContextToID(1,right_context)<<
                                         "] + IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,left_context)<<"]["<<j<<"]["<<j+k<<"]["<< convertContextToID(1,right_context)<<
                                         "] * IT[" <<sizeLeftContext*(*itRightS).id + convertContextToID(0,left_context)<<"]["<<j+k+1<<"]["<<j+i<<"]["<< convertContextToID(1,right_context) << "] * " <<  "ProbR = "
                                         <<p[sizeLeftContext*nonterminal.id+ convertContextToID(0,left_context)][j][i+j][convertContextToID(1,right_context)]<< std::endl;*/
                                //std::cout << p[sizeLeftContext*nonterminal.id+ convertContextToID(0,left_context)][j][i+j][convert_context_to_id(1,right_context)]  << " + " << bInside << " * " << cInside << " * " << (*itRight).second.first << std::endl;
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

Grammar::Grammar::~Grammar() = default;

void Grammar::Grammar::gibbs_sampling_pcfg(int iterations) {
    //srand((unsigned) time(nullptr));
    size_t sentencesLength = words.size();
    for (int j = 0; j < iterations; j++) {
        parse_trees_vec.clear();
        for (size_t i = 0; i< sentencesLength; i++) {
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> vProds;
            double ***insideTable;
            insideTable = cyk_prob_vec(words[i]);
            sample_parse_tree_vec(vProds, rules[0], words[i], insideTable, 0, words[i].size() - 1);
            free_inside_table(insideTable, words[i].size());
            std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> pairTree;
            pairTree.first = words[i];
            pairTree.second = vProds;
            parse_trees_vec.push_back(pairTree);
        }
        std::vector<Symbol::Symbol> a;
        calculate_new_theta_vec_opt(-1);
        update_parse_tress_theta();
        if (j%(iterations/10) == 0) {
            std::pair<double,  double> pMpP = perplexity(words);
            std::cout << "iteration: " << j <<" - PerplexityM: "<< pMpP.first << " - PerplexityP: "<< pMpP.second << std::endl;
        }
    }
}

void Grammar::Grammar::gibbs_sampling_pcsg(int iterations) {
    auto totalTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto insideTTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto iterationTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto newThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto updateThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto perplexityTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    //srand((unsigned) time(nullptr));
    size_t sentencesLength = words.size();
    //auto startMetropolis = std::chrono::system_clock::now();
    for (int j = 0; j < iterations; j++) {
        parse_trees_vec.clear();
        for (size_t i = 0; i< sentencesLength; i++) {
            auto startIt = std::chrono::system_clock::now();
            actual_production.clear();
            actual_production.push_back(non_terminals[0]);
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> vProds;
            auto startInsideTable = std::chrono::system_clock::now();
            double ****insideTableKL;
            insideTableKL = cyk_prob_kl_vec(words[i]);
            insideTTime += std::chrono::system_clock::now() - startInsideTable;
            sample_parse_tree_kl_vec(vProds, rules[0], words[i], insideTableKL, 0, words[i].size() - 1);
            free_inside_table_kl(insideTableKL, words[i].size());
            std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> pairTreeVec;
            pairTreeVec.first = words[i];
            pairTreeVec.second = vProds;
            parse_trees_vec.push_back(pairTreeVec);

            iterationTime = std::chrono::system_clock::now() - startIt;
            //std::cout << "   insideT: " << insideTTime.count() << " iterationT: " << iterationTime.count() << std::endl;
        }
        std::vector<Symbol::Symbol> a;
        auto startNewTheta = std::chrono::system_clock::now();
        calculate_new_theta_vec_opt(-1);
        newThetaTime += std::chrono::system_clock::now() - startNewTheta;
        auto startUpdateTheta = std::chrono::system_clock::now();
        update_parse_tress_theta();
        updateThetaTime += std::chrono::system_clock::now() - startUpdateTheta;

        //std::cout << "   newTheta: " << newThetaTime.count() << " updateTheta: " << updateThetaTime.count() << std::endl;

        if (j%(iterations/5) == 0) {
            //auto startPerplexityTime = std::chrono::system_clock::now();
            //std::pair<double,  double> pMpP = perplexity_kl(words, true);
            //perplexityTime += std::chrono::system_clock::now() - startPerplexityTime;
            //std::cout << "      iteration: " << j <<" - Perplexity: "<< pMpP.first << " - PerplexityN: "<< pMpP.second <<std::endl;
            std::cout << "      iteration: " << j  << std::endl;
        }
    }
}

void Grammar::Grammar::free_inside_table(double ***p, size_t wSize) const {
    for (unsigned long i = 0; i < non_terminals.size(); i++) {
        for (size_t j = 0; j < wSize; j++)
            delete[] p[i][j];
        delete[] p[i];
    }
    delete[] p;
}

void Grammar::Grammar::free_inside_table_kl(double ****p, size_t wSize) const {
    for (size_t i = 0; i < non_terminals.size() * context_amount.first; i++) {
        for (size_t j = 0; j < wSize; j++) {
            for (size_t k = 0; k < wSize; k++)
                delete[] p[i][j][k];
            delete[] p[i][j];
        }
        delete[] p[i];
    }
    delete p;
}

void Grammar::Grammar::sample_parse_tree_kl_vec_opt(
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> &vr,
        Rule::Rule r, std::vector<Symbol::Symbol> w, double ****inside_table, size_t i, size_t k) {
    size_t jRange = k - i;
    size_t rightRange = n_non_terminals * n_non_terminals;

    Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
    std::vector<Symbol::Symbol>::iterator itLeft;
    std::vector<Symbol::Symbol> leftContext;
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
    std::vector<Symbol::Symbol> rightContext;
    while (itLeft != r.left.end() ) {
        rightContext.push_back((*itLeft));
        itLeft++;
    }
    //double sumJBC = 0.0;
    if (jRange > 0) {
        //std::cout << w << std::endl;
        auto *jbc = new double[(jRange+1) * rightRange];

        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight ++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
            for (itSymbol = (*itRight).first.begin(); itSymbol != (*itRight).first.end(); itSymbol++) {
                if (!(*itSymbol).terminal && !(*itSymbol).context) {
                    std::vector<Symbol::Symbol> rightContextLine;
                    itSymbol++;
                    if (context_size.second > 0) {
                        rightContextLine.push_back((*itSymbol));
                        rightContextLine.front().context = true;
                        if (!rightContext.empty()) {
                            rightContextLine.insert(rightContextLine.end(), rightContext.begin(), rightContext.end());
                            rightContext.pop_back();
                        }
                    }
                    itSymbol--;
                    for (unsigned int l = 0; l < jRange; l++) {
                        size_t indexB = (*itSymbol).id * n_non_terminals * jRange;
                        double pBij = inside_table[context_amount.first * (*itSymbol).id +
                                convert_context_to_id(0, leftContext)]
                        [i][i + l][convert_context_to_id(1, rightContextLine)];
                        itSymbol++;
                        double pCjk = inside_table[context_amount.first * (*itSymbol).id +
                                convert_context_to_id(0, leftContext)]
                        [i + l + 1][k][convert_context_to_id(1, rightContext)];
                        size_t indexC = (*itSymbol).id * jRange;
                        jbc[indexB + indexC + l] = ((*itRight).second.first * pBij * pCjk) / inside_table
                        [context_amount.first * nonterminal.id + convert_context_to_id(0, leftContext)]
                        [i][k][convert_context_to_id(1, rightContext)];
                        itSymbol--;
                        /*std::cout << "SumJBC = " <<sumJBC << " ";

                        std::cout << "jbc[" <<indexB+indexC+l<<"] = " << (*itRight).second.first <<
                        " * " << pBij << " ("<< "IT["<<
                                      left_context.size()*(*itSymbol).id+ convertContextToID(0,left_context)<<"]["<<i<<"]["<<i+l<<"]["<<
                                      convertContextToID(1,right_context)<<"])"
                        << " * " << pCjk
                                << " ("<< "IT["<<
                                left_context.size()*(*itSymbol).id+ convertContextToID(0,left_context)<<"]["<<i+l+1<<"]["<<k<<"]["<<
                                convert_context_to_id(1,right_context)<<"])"
                        << " / " << inside_table
                        [left_context.size()*nonterminal.id+ convertContextToID(0,left_context)]
                        [i][k][convertContextToID(1,right_context)] << " IT["<<
                        left_context.size()*nonterminal.id+ convertContextToID(0,left_context)<<"]["<<i<<"]["<<k<<"]["<<
                        convertContextToID(1,right_context)<<"]"<<std::endl;*/


                    }
                    break;
                }
            }
        }


        size_t l = 0;
        for (l = 1; l < rightRange*jRange; l++)
            jbc[l] = jbc[l] + jbc[l-1];
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        //std::cout << "Sensitive: " << jbc[l-1] << std::endl;
        for (l = 0; l < rightRange*jRange; l++)
            if (p < jbc[l])
                break;
        size_t j = 0;

        delete[]jbc;
        Symbol::Symbol rightB = Symbol::Symbol("", 0,false);
        Symbol::Symbol rightC = Symbol::Symbol("", 0,false);

        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        rightProduction = r.get_right_side_by_id(l / (n_non_terminals * jRange),
                                                 (l % (n_non_terminals * jRange)) / jRange, false, n_non_terminals);
        std::vector<Symbol::Symbol>::iterator itSymbol;
        for (itSymbol = rightProduction.first.begin(); itSymbol != rightProduction.first.end(); itSymbol++) {
            if (!(*itSymbol).terminal && !(*itSymbol).context) {
                rightB = (*itSymbol).clone();
                j = l - n_non_terminals * jRange * (*itSymbol).id;
                itSymbol++;
                rightC = (*itSymbol).clone();
                if ((*itSymbol).id > 0)
                    j = j % (*itSymbol).id;
                break;
            }
        }
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);

        vr.push_back(production);
        apply_production(production);
        std::vector<Symbol::Symbol> leftContext2;
        std::vector<Symbol::Symbol> rightContext2;
        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();

        size_t lc , rc;
        std::uniform_int_distribution<size_t> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<size_t> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = ((int) rand()%(leftContext2.size()+1));
        //rc = (int) rand()%(rightContext2.size()+1);
        for (size_t i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (size_t i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::vector<Symbol::Symbol> lhs;
        lhs.insert(lhs.end(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightB);
        lhs.insert(lhs.end(), rightContext2.begin(), rightContext2.end());
        Rule::Rule nextRule = find_rule_by_lhs(lhs);
        sample_parse_tree_kl_vec_opt(vr, nextRule, w, inside_table, i, i + j);

        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();


        lc = distIntL(mt);
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (unsigned int i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (unsigned int i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();

        lhs.clear();
        lhs.insert(lhs.begin(), leftContext2.begin(), leftContext2.end());
        lhs.push_back(rightC);
        lhs.insert(lhs.begin(), rightContext2.begin(), rightContext2.end());
        sample_parse_tree_kl_vec_opt(vr, find_rule_by_lhs(lhs), w, inside_table, i + j + 1, k);

    } else {
        std::vector<Symbol::Symbol> leftContext2;
        std::vector<Symbol::Symbol> rightContext2;
        get_actual_context(leftContext2, rightContext2);

        while (leftContext2.size() > context_size.first)
            leftContext2.erase(leftContext2.begin());
        while (rightContext2.size() > context_size.second)
            rightContext2.pop_back();
        size_t lc , rc;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<size_t> distIntL(0, leftContext2.size());
        lc = distIntL(mt);
        std::uniform_int_distribution<size_t> distIntR(0, rightContext2.size());
        rc = distIntR(mt);
        //lc = (int) rand()%(leftContext2.size()+1);
        //rc = (int) rand()%(rightContext2.size()+1);
        for (size_t i2 = 0; i2 < lc; i2++)
            leftContext2.erase(leftContext2.begin());
        for (size_t i2 = 0; i2 < rc; i2++)
            rightContext2.pop_back();
        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightProduction;
        rightProduction.first.insert(rightProduction.first.end(), leftContext.begin(), leftContext.end());
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        std::string s = w[i].name;
        for (itRight = r.right.begin(); itRight != r.right.end(); itRight++) {
            std::vector<Symbol::Symbol>::iterator itSymbol;
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
        std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> production = std::make_pair(r.left, rightProduction);
        vr.push_back(production);
        apply_production(production);
        return;
    }
}

/*
double ****Grammar::CYKProbKLVecOpt(std::vector<Symbol::Symbol> w) {
    int sizeLeftContext = 0;
    int sizeRightContext = 0;
    for (unsigned long i = 0; i <= context_size.first; i++) {
        std::vector<Symbol::Symbol> word;
        std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
        generatePermutation(permutationsTerminals, terminals, i, word, true);
        sizeLeftContext += permutationsTerminals.size();

    }
    context_amount.first = sizeLeftContext;
    for (unsigned long i = 0; i <= context_size.second; i++) {
        std::vector<Symbol::Symbol> word;
        std::vector<std::vector<Symbol::Symbol>> permutationsTerminals;
        generate_permutation(permutationsTerminals, non_terminals, i, word, true);
//        rightContexts.insert(rightContexts.end(), permutationsTerminals.begin(), permutationsTerminals.end());
        sizeRightContext += permutationsTerminals.size();
    }
    context_amount.second = sizeRightContext;
    auto ****p = new double***[non_terminals.size()*sizeLeftContext];
    for (unsigned long i = 0; i < non_terminals.size()*sizeLeftContext; i++) {
        p[i] = new double **[w.size()];
        for (unsigned long j = 0; j < w.size(); j++) {
            p[i][j] = new double *[w.size()];
            for (unsigned long k = 0; k < w.size(); k++)
                p[i][j][k] = new double[sizeRightContext];
        }
    }

    for (unsigned long i = 0; i < non_terminals.size()*sizeLeftContext; i++)
        for (unsigned long j = 0; j < w.size(); j++)
            for (unsigned long k = 0; k < w.size(); k++)
                for (int l = 0; l < sizeRightContext; l++)
                    p[i][j][k][l] = 0.0;
    std::vector<Rule>::iterator itRule;
    for (unsigned long i = 0; i < w.size(); i++) {
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<Symbol::Symbol>::iterator itLeft;
            Symbol::Symbol nonterminal =  Symbol::Symbol("", 0,false);
            std::vector<Symbol::Symbol> left_context;
            for (itLeft = (*itRule).left.begin(); itLeft != (*itRule).left.end(); itLeft++) {
                if (!itLeft->terminal && !itLeft->context) {
                    nonterminal = (*itLeft).clone();
                    itLeft++;
                    break;
                }
                else {
                    left_context.push_back((*itLeft));
                }
            }
            std::vector<Symbol::Symbol> right_context;
            while (itLeft != (*itRule).left.end() ) {
                right_context.push_back((*itLeft));
                itLeft++;
            }

            std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightHandSide;
            rightHandSide = (*itRule).get_right_side_by_id(w[i].id, 0, true, n_non_terminals);
            std::vector<Symbol::Symbol>::iterator itRightS;
            for (itRightS = rightHandSide.first.begin(); itRightS != rightHandSide.first.end(); itRightS++) {
                if (itRightS->terminal && !itRightS->context) {
                    if(itRightS->name == w[i].name) {
                        p[sizeLeftContext*nonterminal.id+ convert_context_to_id(0,left_context)][i][i][convertContextToID(1,right_context)] = rightHandSide.second.first;
                        */
/*std::cout<< "IT["<<sizeLeftContext*nonterminal.id + convert_context_to_id(0,left_context)<<"]["<<i<<"]["<<i<<"]["<< convertContextToID(1,right_context)<<"] = ProbR"
                        << " " << (*itRight).second.first << std::endl;*//*

                        break;
                    }
                }
            }

        }
    }
    for (unsigned long i = 1; i < w.size(); i++) {
        for (unsigned long j = 0; j < w.size()-i; j++) {
            for (unsigned long k = 0; k< i; k++) {
                for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
                    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
                    for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                            if (!(*itRight).first[(*itRule).index_1st_non_context].terminal &&
                                !(*itRight).first[(*itRule).index_1st_non_context].context) {
                                std::vector<Symbol::Symbol> rightContextLine;
                                if ( context_size.second > 0) {
                                    rightContextLine.push_back((*itRight).first[(*itRule).index_1st_non_context+1]);
                                    rightContextLine.front().context = true;
                                    if (!(*itRule).right_context.empty()) {
                                        rightContextLine.insert(rightContextLine.end(), (*itRule).right_context.begin(), (*itRule).right_context.end());
                                        rightContextLine.pop_back();
                                    }
                                }
                                int leftContextID = convertContextToID(0,(*itRule).left_context);
                                int rightContextID = convert_context_to_id(1,(*itRule).right_context);
                                double bInside = p[sizeLeftContext*(*itRight).first[(*itRule).index_1st_non_context].id + leftContextID][j][j+k][convertContextToID(1,rightContextLine)];
                                double cInside = p[sizeLeftContext*(*itRight).first[(*itRule).index_1st_non_context+1].id + leftContextID][j+k+1][j+i][rightContextID];

                                p[sizeLeftContext*(*itRight).first[(*itRule).index_1st_non_context].id+ leftContextID][j][i+j][rightContextID] =
                                        p[sizeLeftContext*(*itRight).first[(*itRule).index_1st_non_context].id+ leftContextID][j][i+j][rightContextID] + bInside*cInside*(*itRight).second.first;



*/
/*                                for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
                                    std::vector<Symbol::Symbol>::iterator itRightS;
                                    std::vector<Symbol::Symbol> rightContextLine;
                                    if ( context_size.second > 0) {
                                        rightContextLine.push_back((*itRight).first[(*itRule).index_1st_non_context+1]);
                                        rightContextLine.front().context = true;
                                        if (!(*itRule).right_context.empty()) {
                                            rightContextLine.insert(rightContextLine.end(), (*itRule).right_context.begin(), (*itRule).right_context.end());
                                            rightContextLine.pop_back();
                                        }
                                    }
                                    int leftContextID = convertContextToID(0,(*itRule).left_context);
                                    int rightContextID = convertContextToID(1,(*itRule).right_context);
                                    double bInside = p[sizeLeftContext*((*itRight).first[(*itRule).index_1st_non_context]).id + leftContextID][j][j+k][convert_context_to_id(1,rightContextLine)];
                                    double cInside = p[sizeLeftContext*((*itRight).first[(*itRule).index_1st_non_context+1]).id + leftContextID][j+k+1][j+i][rightContextID];
                                    //nonterminal  é igual a (*itRight).first[(*itRule).index_1st_non_context]
                                    p[sizeLeftContext*(*itRight).first[(*itRule).index_1st_non_context].id+ leftContextID][j][i+j][rightContextID] =
                                            p[sizeLeftContext*(*itRight).first[(*itRule).index_1st_non_context].id+ leftContextID][j][i+j][rightContextID] + bInside*cInside*(*itRight).second.first;
                                    break;
                                }*//*


                        }
                    }
                }
            }
        }
    }
    //printInsideTableKL(p,w.size());
    return p;
}
*/


void Grammar::Grammar::generate_rules_regular() {
    for (const auto& n: non_terminals) {
        std::vector<Symbol::Symbol> left;
        left.push_back(n);
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right;
        for (const auto& t: terminals) {
            std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rightHS;
            rightHS.first.push_back(t);
            for (const auto& nr: non_terminals) {
                rightHS.first.push_back(nr);
                rightHS.second.first = 1.0/ static_cast<double>(n_non_terminals * n_terminals + 1);
                rightHS.second.second = ALFA;
                right.push_back(rightHS);
                rightHS.first.pop_back();
            }
        }
        std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rightHS;
        rightHS.second.first = 1.0/ static_cast<double>(n_non_terminals * n_terminals + 1);
        rightHS.second.second = ALFA;
        right.push_back(rightHS);
        Rule::Rule r = Rule::Rule(left, right);
        rules.push_back(r);
    }
}

void Grammar::Grammar::baum_welch(int iteration) {
    std::vector<double> countNd;
    std::vector<double> countNl;
    vector<vector<vector<double>>> countNT;
    for (int i = 0; i < iteration; i++) {
        cout << "Iteration: "<<i << endl;
        baum_welch_expectation(countNd, countNl, countNT);
        baum_welch_maximization(countNd, countNl, countNT);
    }

}


void Grammar::Grammar::baum_welch_expectation(std::vector<double> &count_nd, std::vector<double> &count_nl, std::vector<std::vector<std::vector<double>>> &count_nt) {
    count_nd.clear();
    count_nl.clear();
    count_nt.clear();
    for (const auto& n: non_terminals) {
        count_nd.push_back(0);
        count_nl.push_back(0);
        count_nt.emplace_back();
        for (const auto& t: terminals) {
            count_nt[n.id].push_back(vector<double>());
            for (auto nl: non_terminals) {
                count_nt[n.id][t.id].push_back(0.0);
            }
        }
    }
    for (auto x: words) {
        vector<vector<double>> F, B;
        for (unsigned long i = 0; i <= x.size(); i++) {
            F.emplace_back();
            B.emplace_back();
            for (auto n: non_terminals) {
                F[i].push_back(0.0);
                B[i].push_back(0.0);
            }
        }
        //To DO: contar quantas vezes x acontece ou noã precisa?
        for (const auto& n: non_terminals) {
            if (n.id == 0) F[0][n.id] = 1; else F[0][n.id] = 0;
            B[x.size()][n.id] = rules[n.id].right[n_non_terminals * n_terminals].second.first;
            //cout << "F[" << "0" << "][" << n.id << "] = " << F[0][n.id] << "\t\t\t\t\t\t\t\t";
            //cout << "B[" << x.size() << "][" << n.id << "] = "
            //     << rules[n.id].right[n_non_terminals * n_terminals].second.first << endl;
        }
        for (unsigned long i = 1; i<= x.size(); i++) {
            for (const auto& n: non_terminals) {
                F[i][n.id] = 0;
                B[x.size()-i][n.id]=0;
                for (const auto& nl: non_terminals) {
                    F[i][n.id] += F[i-1][nl.id]*rules[nl.id].right[x[i-1].id * n_terminals + n.id].second.first;
                //    cout << "F["<<i<<"]"<<"["<<n.id<<"]"<< " += F["<<i-1<<"]"<<"["<<nl.id<<"] * " << "P["<<nl.id<<"]["<<x[i-1].id<<"]["<<n.id<<"] = "
                //    << F[i-1][nl.id] << "*" << rules[nl.id].right[x[i-1].id*n_terminals+n.id].second.first << " = " << F[i][n.id]<< "\t\t\t\t";
                    B[x.size()-i][n.id] +=  B[x.size()-i+1][nl.id]* rules[n.id].right[x[x.size()-i].id * n_terminals + nl.id].second.first;
                //    cout << "B["<<x.size()-i<<"]"<<"["<<n.id<<"]"<< " += " << "B["<<x.size()-i+1<<"]"<<"["<<nl.id<<"]" << " * P["<<n.id<<"]["<<x[i-1].id<<"]["<<nl.id<<"] = "
                //    << B[x.size()-i+1][nl.id]<< "*" << rules[n.id].right[x[x.size()-i].id*n_terminals+nl.id].second.first << " = " << B[x.size()-i][n.id] <<endl;
                }
            }
        }
        double total = 0.0;
        //double totalb = B[0][0];
        for (const auto& n: non_terminals) {
            total += F[x.size()][n.id]*rules[n.id].right[n_non_terminals * n_terminals].second.first;
            //cout << " total = " << "total + " <<  F[x.size()][n.id] << "*" << rules[n.id].right[n_non_terminals*n_terminals].second.first << endl;
        }
        //cout << "Total = " << total << endl;
        //cout << "TotalB = " << totalb <<endl;
        for (unsigned long i = 1; i<= x.size(); i++) {
            for (const auto& n: non_terminals) {
                for (const auto& nl: non_terminals) {
                    double val = (F[i-1][n.id] * rules[n.id].right[x[i-1].id * n_terminals + nl.id].second.first * B[i][nl.id]) / total;
                    //cout << "val = F[" <<i-1<<"]["<<n.id<<"]* P["<<n.id<<"]["<<x[i-1].id<<"]["<<nl.id<<"]*"<<"B["<<i<<"]["<<nl.id<<"] = " <<
                    //F[0][0] << "*" << rules[n.id].right[x[i-1].id*n_terminals+nl.id].second.first <<"*" << B[i][nl.id]<< " = " << val << " / ";
                    count_nt[n.id][x[i - 1].id][nl.id] += val;
                    //cout << "count_nt["<<n.id<<"]["<<x[i-1].id<<"]["<<nl.id<<"] += val  = " <<count_nt[n.id][x[i-1].id][nl.id] << " / ";
                    count_nl[nl.id] += val;
                    //cout << "count_nl["<<nl.id<<"] += "<< count_nl[nl.id] << endl;
                }

            }
        }
        for (const auto& n: non_terminals) {
            count_nd[n.id] += (F[x.size()][n.id] * rules[n.id].right[n_non_terminals * n_terminals].second.first) / total;
            //cout << "count_nd["<<n.id<<"] += F["<<x.size()<<"]["<<n.id<<"]"<<"*"<< "P["<<n.id<<"][F]["<<n_non_terminals*n_terminals<<"] = "
            //<< F[x.size()][n.id] << "*" <<rules[n.id].right[n_non_terminals*n_terminals].second.first << " / " << total << " = " << count_nd[n.id] << endl;
            if (n.id == 0) count_nl[n.id] += B[0][n.id] / total;
        }
    }
}

void Grammar::Grammar::baum_welch_maximization(vector<double> &count_nd_r, vector<double> &count_nl_r,
                                      vector<std::vector<std::vector<double>>> &count_ntr) {
    vector<double> val;
    for (auto n: non_terminals) {
        val.push_back(0.0);
    }
    for (const auto& n: non_terminals) {
        rules[n.id].right[n_non_terminals * n_terminals].second.first = (1.0 * count_nd_r[n.id]) / (1.0 * count_nl_r[n.id]);
        //cout << "rulesF["<<n.id<<"]["<<n_non_terminals*n_terminals<<"] = count_nd_r["<<n.id<<"]/count_nl_r["<<n.id<<"] = " << count_nd_r[n.id] <<"/"<< count_nl_r[n.id] << endl;
        double total = 0.0;
        for (const auto& t: terminals) {
            for (const auto& nl: non_terminals) {
                total += 1.0 * count_ntr[n.id][t.id][nl.id];
                rules[n.id].right[t.id * n_terminals + nl.id].second.first = (1.0 * count_ntr[n.id][t.id][nl.id]) / (1.0 * count_nl_r[n.id]);
                //cout << "rulesF["<<n.id<<"]["<<t.id*n_terminals+nl.id<<"] = count_ntr["<<n.id<<"]["<<t.id<<"]["<<nl.id<<"]/count_nl_r["<<n.id<<"] = "
                //<< count_ntr[n.id][t.id][nl.id] <<"/"<< count_nl_r[n.id] << endl;
            }
        }
        //cout << "cNdR: " << count_nd_r[n.id] << " cNTR: " << total
        //<< " sum: " << count_nd_r[n.id]+ total<< " cNlR: "<< count_nl_r[n.id] << endl;
    }
}

void Grammar::Grammar::inside_outside(int iterations) {
    for (int i = 0; i < iterations; i++) {
        vector<vector<vector<double>>>newProbsNTn;
        vector<vector<vector<double>>>newProbsNTd;
        vector<vector<double>>newProbsTn;
        for (size_t i2 = 0 ; i2 < n_non_terminals; i2++) {
            newProbsNTn.emplace_back();
            newProbsNTd.emplace_back();
            for (size_t j = 0; j < n_non_terminals; j++) {
                newProbsNTn[i2].push_back(vector<double>());
                newProbsNTd[i2].push_back(vector<double>());
                for (size_t k = 0; k < n_non_terminals; k++) {
                    newProbsNTn[i2][j].push_back(0.0);
                    newProbsNTd[i2][j].push_back(0.0);
                }
            }
        }
        for (size_t i2 = 0 ; i2 < n_non_terminals; i2++) {
            newProbsTn.emplace_back();
            for (size_t m = 0; m < n_terminals; m++) {
                newProbsTn[i2].push_back(0.0);
            }
        }

        cout << "Iteration: "<<i << endl;
        for (auto w: words) {
            double *** insideTable;
            insideTable = cyk_prob_vec(w);
            //print_inside_table(insideTable, w.size());
            std::vector<std::vector<std::vector<double>>> oT = outside_table(w, insideTable);
            //print_outside_table(oT);
            vector<vector<vector<vector<vector<double>>>>> W;
            vector<vector<vector<double>>> V;
            calculate_wv(W, V, w, insideTable, oT);

            for (unsigned long s = 0 ; s <w.size()-1; s++) {
                for (unsigned long t = s+1; t <w.size(); t++) {
                    for (size_t i2 = 0; i2 < n_non_terminals; i2++) {
                        for (size_t j = 0; j < n_non_terminals; j++) {
                            for (size_t k = 0; k < n_non_terminals; k++) {
                                    newProbsNTn[i2][j][k] += W[s][t][i2][j][k];
                                    //cout << "nPNTn["<<i2<<"]["<<j<<"]["<<k<<"] = "<< newProbsNTn[i2][j][k] << endl;
                            }
                        }
                    }
                }
            }

            for (size_t i2 = 0; i2 < n_non_terminals; i2++) {
                for (size_t j = 0; j < n_non_terminals; j++) {
                    for (size_t k = 0; k < n_non_terminals; k++) {
                        for (unsigned long s = 0 ; s <w.size(); s++) {
                            for (unsigned long t = s; t <w.size(); t++) {
                                    newProbsNTd[i2][j][k] += V[s][t][i2];
                                    //cout << "nPNTd["<<i2<<"]["<<j<<"]["<<k<<"] = "<< newProbsNTd[i2][j][k] << endl;
                            }
                        }
                    }
                }
            }

            for (unsigned long t = 0; t <w.size(); t++) {
                for (size_t i2 = 0; i2 < n_non_terminals; i2++) {
                    for (size_t m = 0; m < n_terminals; m++) {
                        if(terminals[m].id == w[t].id) {
                            newProbsTn[i2][m] += V[t][t][i2];
                            //cout << "nPT["<<i2<<"]["<<m<<"] = "<< newProbsTn[i2][m] << endl;
                        }
                    }
                }
            }
        }

        vector<Rule::Rule>::iterator itRule;
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
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
        print_grammar();
    }
}

void Grammar::Grammar::calculate_wv(vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &w,
                           vector<std::vector<std::vector<double>>> &v, const vector<Symbol::Symbol>& word, double***inside_table, vector<vector<vector<double>>>o_t) {

    for (unsigned long s = 0 ; s < word.size(); s++) {
        w.emplace_back();
        for (unsigned long t = 0; t < word.size(); t++) {
            w[s].push_back(vector<vector<vector<double>>>());
            for (size_t i = 0; i < n_non_terminals; i++) {
                w[s][t].push_back(vector<vector<double>>());
                for (size_t j = 0; j < n_non_terminals; j++) {
                    w[s][t][i].push_back(vector<double>());
                    for (size_t k = 0; k < n_non_terminals; k++) {
                        w[s][t][i][j].push_back(0.0);
                    }
                }
            }
        }
    }

    for (unsigned long s = 0 ; s < word.size(); s++) {
        v.emplace_back();
        for (unsigned long t = 0; t < word.size(); t++) {
            v[s].push_back(vector<double>());
            for (size_t i = 0; i < n_non_terminals; i++) {
                v[s][t].push_back(0.0);
            }
        }
    }

    for (unsigned long s = 0 ; s < word.size(); s++) {
        for (unsigned long t = 0; t < word.size(); t++) {
            for (size_t i = 0; i < n_non_terminals; i++) {
                for (size_t j = 0; j < n_non_terminals; j++) {
                    for (size_t k = 0; k < n_non_terminals; k++) {
                        for (unsigned long r = s; r < t; r++) {
                            w[s][t][i][j][k] += (rules[i].right[k + n_non_terminals * j].second.first *
                                                 inside_table[j][s][r] * inside_table[k][r + 1][t] * o_t[i][s][t]) / inside_table[0][0][word.size() - 1];
                            //cout<< "w["<<s<<"]["<<t<<"]["<<i<<"]["<<j<<"]["<<k<<"] = " << w[s][t][i][j][k] << endl;
                        }
                    }
                }
            }
        }
    }

    for (unsigned long s = 0 ; s < word.size(); s++) {
        for (unsigned long t = 0; t < word.size(); t++) {
            for (size_t i = 0; i < n_non_terminals; i++) {
                v[s][t][i] = (inside_table[i][s][t] * o_t[i][s][t]) / inside_table[0][0][word.size() - 1];
                //cout<< "v["<<s<<"]["<<t<<"]["<<i<<"] = " << v[s][t][i] << endl;
            }
        }
    }
}

void Grammar::Grammar::gen_fpta() {
    Symbol::Symbol nI = Symbol::Symbol("NTI", 0, false, false);
    non_terminals.push_back(nI);
    int nNT = 1;
    vector<Symbol::Symbol> vr;
    vr.push_back(nI);
    Rule::Rule r = Rule::Rule(vr, std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>());
    rules.push_back(r);
    for (const auto& w: words) {
        string s;
        for (auto & i : w) {
            s += i.name;
            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
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
                Symbol::Symbol nt = Symbol::Symbol("NT"+s, nNT, false, false);
                non_terminals.push_back(nt);
                nNT++;
                std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> right;
                right.first.push_back(i);
                right.first.push_back(nt);
                right.second.first = 1.0;
                rules[nI.id].right.insert(rules[nI.id].right.end()-1, right);
                //rules[nI.id].right.push_back(right);
                vr.clear();
                vr.push_back(nt);
                r = Rule::Rule(vr, std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>());
                rules.push_back(r);
                nI = nt;
            }

        }
        bool existRule = false;
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = rules[nI.id].right.begin(); itRight != rules[nI.id].right.end(); itRight++) {
            if ((*itRight).first[0].name.empty()) {
                (*itRight).second.first += 1.0;
                existRule = true;
            }
        }
        if (!existRule) {
            Symbol::Symbol nt = Symbol::Symbol("", -1, false, false);
            non_terminals.push_back(nt);
            std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> right;
            right.first.push_back(nt);
            right.second.first = 1.0;
            rules[nI.id].right.push_back(right);

        }
        nI = non_terminals[0];
    }
    n_non_terminals = nNT;
}

void Grammar::Grammar::alergia(double alpha) {
    non_terminals.clear();
    rules.clear();
    n_non_terminals = 0;
    gen_fpta();
    print_grammar();
    vector<Symbol::Symbol> red;
    red.push_back(non_terminals[0]);
    vector<Symbol::Symbol> blue;
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
            if (compatible_alergia(cRed, (*itBlue), alpha)) {
                cout << "merge " << cRed.name << " and " << (*itBlue).name << endl;
                stochastic_merge(cRed, (*itBlue));
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
        print_grammar();

    }
    remove_unused_nt();
    normalize_probs();
}

bool Grammar::Grammar::compatible_alergia(const Symbol::Symbol& a, const Symbol::Symbol& b, double alpha) {
    bool correct = true;
    if (!test_alergia((*(rules[a.id].right.end() - 1)).second.first, rules[a.id].freq(),
                      (*(rules[b.id].right.end() - 1)).second.first, rules[b.id].freq(), alpha))
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
        if (!test_alergia(dFreqa, rules[a.id].freq(), dFreqb, rules[b.id].freq(), alpha))
            correct = false;
    }
    return correct;
}

bool Grammar::Grammar::test_alergia(double probFa, double freqa, double probFb, double freqb, double alpha) {
    double y = abs(probFa/freqa - probFb/freqb);
    return (y < ( (sqrt(1/freqa) + sqrt(1/freqb)) * sqrt(0.5 * log(2/alpha)) ));
}

void Grammar::Grammar::stochastic_merge(const Symbol::Symbol& a, const Symbol::Symbol& b) {
    vector<Rule::Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule < rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            if ((*itRight).first[1].id == b.id) {
                stochastic_fold(a, b);
                std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> newRight;
                //int n = (*itRight).second.first;
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

void Grammar::Grammar::stochastic_fold(const Symbol::Symbol& a, const Symbol::Symbol& b) {
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = rules[b.id].right.begin(); itRight != rules[b.id].right.end()-1; itRight++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRighta;
        for (itRighta = rules[a.id].right.begin(); itRighta != rules[a.id].right.end(); itRighta++) {
            if (itRighta->first[0].id == itRight->first[0].id) {
                if (itRighta->second.first >=1.0)
                    break;
            }
        }
        if (itRighta != rules[a.id].right.end()) {
            if (itRighta->second.first >= 1.0) {
                stochastic_fold(itRighta->first[1], itRight->first[1]);
                (*itRighta).second.first += (*itRight).second.first;
            }
        }
        else {
            itRighta->first[1] = Symbol::Symbol(itRight->first[1]);
            (*itRighta).second.first = (*itRight).second.first;

        }
    }
    (*(rules[a.id].right.end()-1)).second.first += (*(rules[b.id].right.end()-1)).second.first;

}

void Grammar::Grammar::remove_unused_nt() {
    vector<Symbol::Symbol> notUsed;
    vector<Symbol::Symbol> used;
    used.push_back(rules[0].left[0]);
    for (unsigned long i = 0; i < used.size(); i++) {
        Rule::Rule r = rules[used[i].id];
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = r.right.begin(); itRight != r.right.end()-1; itRight++) {
            if ((*itRight).second.first <1.0) {
                notUsed.push_back((*itRight).first[1]);
                Symbol::Symbol aux = (*itRight).first[1];
                recursive_insert_unused(notUsed, aux);

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
        vector<Rule::Rule>::iterator itRule;
        for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
            if ((*itRule).left[0].id == nu.id) {
                rules.erase(itRule);
                break;
            }
        }
    }
}

void Grammar::Grammar::recursive_insert_unused(vector<Symbol::Symbol> &unused, const Symbol::Symbol& nt) {
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = rules[nt.id].right.begin();itRight != rules[nt.id].right.end()-1; itRight++) {
        unused.push_back((*itRight).first[1]);
        recursive_insert_unused(unused, (*itRight).first[1]);
    }
}

void Grammar::Grammar::normalize_probs() {
    vector<Rule::Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        double freq = itRule->freq();
        for(itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++)
            itRight->second.first = itRight->second.first/freq;
    }

}

void Grammar::Grammar::sample_regular_rules(
        vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>> &vr,
        const std::vector<Symbol::Symbol>& w) {
    Rule::Rule actualState = rules[0];
    pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> prod;
    for (const auto& s: w) {
        double total = 0.0;
        vector<double> probsTerminal;
        for (size_t i = n_terminals * s.id; i < n_terminals * s.id + n_non_terminals; i++) {
            probsTerminal.push_back(total + actualState.right[i].second.first);
            total += actualState.right[i].second.first;
        }
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double p = dist(mt);
        unsigned long i;
        for (i = 0; i < probsTerminal.size(); i++) {
            if (p * total < probsTerminal[i])
                break;
        }

        prod = make_pair(actualState.left, actualState.right[n_terminals * s.id + i]);
        vr.push_back(prod);
        actualState = rules[actualState.right[i].first[1].id];
    }
    prod = make_pair(actualState.left, actualState.right[actualState.right.size()-1]);
    vr.push_back(prod);
}

void Grammar::Grammar::collapsed_gibbs_sample_pfa(int iterations) {
    //srand((unsigned) time(nullptr));
    size_t sentencesLength = words.size();
    for (int j = 0; j < iterations; j++) {
        parse_trees_vec.clear();
        for (size_t i = 0; i< sentencesLength; i++) {
            std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> vProds;
            sample_regular_rules(vProds, words[i]);
            std::pair<std::vector<Symbol::Symbol>, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>>> pairTree;
            pairTree.first = words[i];
            pairTree.second = vProds;
            parse_trees_vec.push_back(pairTree);
        }
        std::vector<Symbol::Symbol> a;
        calculate_new_theta_vec(a);
    }
}

void Grammar::Grammar::generate_n_gram_rules() {
    size_t lastSmallRuleI = 0;
    for (size_t i = 1; i < n_non_terminals; i++) {
        lastSmallRuleI += (int) pow(terminals.size(), i);
    }

    for (size_t i = 0; i <= lastSmallRuleI; i++) {
        vector<Symbol::Symbol> left;
        left.push_back(non_terminals[i]);
        vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
        for (const auto& t: terminals) {
            std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rightHS;
            rightHS.first.push_back(t);
            rightHS.first.push_back(non_terminals[n_terminals * i + t.id + 1]);
            rightHS.second.first = 0;
            rightHS.second.second = ALFA;
            right.push_back(rightHS);
        }
        std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rightHS;
        //rightHS.second.first = 1.0/((n_terminals+1)*1.0);
        rightHS.second.first = 0;
        rightHS.second.second = ALFA;
        right.push_back(rightHS);
        Rule::Rule r = Rule::Rule(left, right);
        rules.push_back(r);
    }

    for (size_t i = lastSmallRuleI+1; i < non_terminals.size(); i++) {
        vector<Symbol::Symbol> left;
        left.push_back(non_terminals[i]);
        vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
        for (const auto& t: terminals) {
            std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rightHS;
            rightHS.first.push_back(t);
            rightHS.first.push_back(non_terminals[lastSmallRuleI + 1 + (((i - lastSmallRuleI - 1) % ((non_terminals.size() - lastSmallRuleI - 1) / n_terminals)) * n_terminals + t.id)]);
            //rightHS.second.first = 1.0/((n_terminals+1)*1.0);
            rightHS.second.first = 0;
            rightHS.second.second = ALFA;
            right.push_back(rightHS);
        }
        std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rightHS;
        rightHS.second.first = 0;
        rightHS.second.second = ALFA;
        right.push_back(rightHS);
        Rule::Rule r = Rule::Rule(left, right);
        rules.push_back(r);
    }

}


void Grammar::Grammar::generate_n_gram_non_terminals() {
    int nIds = 0;
    Symbol::Symbol nt = Symbol::Symbol("NTI", nIds, false, false);
    non_terminals.push_back(nt);
    nt.name = "NT";
    for (size_t i = 1; i <= n_non_terminals; i++)
        recursive_add_terminal_to_nt(nt, i, nIds);
}

void Grammar::Grammar::recursive_add_terminal_to_nt(Symbol::Symbol nt, size_t n, int &nIds) {
    if (n ==0) {
        non_terminals.push_back(nt);
    } else {
        for (const auto& t: terminals) {
            nt.name += t.name;
            nt.id = ++nIds;
            recursive_add_terminal_to_nt(nt, n - 1, nIds);
            nt.name = nt.name.substr(0, nt.name.size()-1);
        }
    }
}

void Grammar::Grammar::train_n_gram() {
    for (auto w: words) {
        vector<Symbol::Symbol> nGram;
        Symbol::Symbol finalLHS = Symbol::Symbol("NTI", 0, false, false);
        for(unsigned long i = 0; i < w.size(); i++) {
            if (i == (unsigned long) n_non_terminals ) break;
            Symbol::Symbol nt = Symbol::Symbol("NTI", 0, false, false);
            if (i >=1)
                nt.name = "NT";
            finalLHS.name = "NT";
            for (const auto& s: nGram) {
                nt.name += s.name;
                finalLHS.name += s.name;
            }
            add_n_gram_rule_frequency(nt, w[i]);
            finalLHS.name += w[i].name;
            nGram.push_back(w[i]);
        }
        for(size_t i = n_non_terminals; i < w.size(); i++) {
            Symbol::Symbol nt = Symbol::Symbol("NT", 0, false, false);
            for (const auto& s: nGram)
                nt.name += s.name;
            add_n_gram_rule_frequency(nt, w[i]);
            finalLHS.name = nt.name;
            nGram.erase(nGram.begin());
            nGram.push_back(w[i]);
        }
        Symbol::Symbol final = Symbol::Symbol("", 0, true, false);
        add_n_gram_rule_frequency(finalLHS, final);
    }
    normalize_probs();
}

void Grammar::Grammar::add_n_gram_rule_frequency(const Symbol::Symbol& lhs, const Symbol::Symbol& next_symbol) {
    vector<Rule::Rule>::iterator itRule;

    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        if (itRule->left[0].name == lhs.name) {
            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
            bool flag = false;
            for (itRight = itRule->right.begin(); itRight != itRule->right.end()-1; itRight++){
                if (itRight->first[0].name == next_symbol.name) {
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
void Grammar::Grammar::prob_sequitur() {
    non_terminals.clear();
    rules.clear();
    n_non_terminals = 0;
    n_non_terminals = 1;
    Symbol::Symbol start = Symbol::Symbol("NT"+ to_string(n_non_terminals-1), n_non_terminals-1, false, false);
    non_terminals.push_back(start);
    vector<Symbol::Symbol> lhs;
    lhs.push_back(start);
    std::unordered_map<string, int> digram_map;
    std::unordered_map<string, tuple<int, int, int>> digram_position;
    std::unordered_map<int, int> rule_map;
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> rhss;
    rules.push_back(Rule::Rule(lhs,rhss));
    for (int i = 0; i < words.size(); i++) {
        //print_grammar();
        cout << "word: " << i <<": " ;
        for (auto a: words[i]) cout << a.name << " ";
        cout << endl;
        if (i == 1049)
            cout <<"go fishing" << endl;

        auto w = words[i];
        std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;
        rhs.first.push_back(Symbol::Symbol(w[0]));
        rhs.second.first = 1.0;
        rhs.second.second = ALFA;
        rhss.push_back(rhs);
        rules[0].right.push_back(rhs);
        int dicard = 0;
        for (auto s: w) {
            //print_grammar();
            if (!(dicard <1)) {
                rules[0].right[i].first.push_back(s);
                string digram = rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " +rules[0].right[i].first[rules[0].right[i].first.size() - 1].name;
                string emptySaux = "";
                string saux = "NT692 NT692";
                string saux2 = "NT1226 NT1226";
                if (saux.compare(digram) == 0)
                    cout<< "Found Digram saux." << endl;
                if (saux2.compare(digram) == 0)
                    cout<< "Found Digram saux2." << endl;
                while (verify_duplicate_digram(digram_position, digram)) {
                    if (saux.compare(digram) == 0)
                        cout<< "Found Digram saux." << endl;
                    if (saux2.compare(digram) == 0)
                        cout<< "Found Digram saux2." << endl;
                    if (digram_position.find(emptySaux) != digram_position.end()) {
                        cout <<" Tem digrama vazio" << endl;
                    }
                    if (i ==3527 || i == 3529) {
                        for (auto m: digram_position) {
                            if (get<2>(m.second) ==3527)
                                cout << "Display m: "<< m.first << ": " << get<1>(m.second) << endl;
                        }
                        if (digram_position.find(saux) != digram_position.end()) {
                            cout << "Achou";
                        }
                        cout << digram << ": " << get<1>(digram_position[digram]) << " : " << endl;
                    }
                    if (!((get<0>(digram_position[digram]) != 0) && (rules[get<0>(digram_position[digram])].right[get<2>(digram_position[digram])].first.size() == 2 ))) {
                        if (rules[0].right[i].first.size() - 1 - get<1>(digram_position[digram]) > 1 || get<2>(digram_position[digram]) != i) {
                            if (!rules[0].right[i].first[rules[0].right[i].first.size() - 2].terminal)
                                rules[rules[0].right[i].first[rules[0].right[i].first.size() - 2].id].right[0].second.first -= 1.0;
                            if (!rules[0].right[i].first[rules[0].right[i].first.size() - 1].terminal)
                                rules[rules[0].right[i].first[rules[0].right[i].first.size() - 1].id].right[0].second.first -= 1.0;
                            enforce_digram_uniqueness(digram_map,  digram_position, i);
                        }
                        else
                            break;
                    } else {
                        string digram_to_look = rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " + rules[0].right[i].first[rules[0].right[i].first.size() - 1].name;
                        if (!rules[0].right[i].first[rules[0].right[i].first.size() - 2].terminal)
                            rules[rules[0].right[i].first[rules[0].right[i].first.size() - 2].id].right[0].second.first -= 1.0;
                        if (!rules[0].right[i].first[rules[0].right[i].first.size() - 1].terminal)
                            rules[rules[0].right[i].first[rules[0].right[i].first.size() - 1].id].right[0].second.first -= 1.0;
                        rules[0].right[i].first.pop_back();
                        //digram_map.erase(rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " + rules[0].right[i].first[rules[0].right[i].first.size() - 1].name);
                        if (rules[0].right[i].first.size() > 1) {
                            string digram_to_erase = rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " + rules[0].right[i].first[rules[0].right[i].first.size() - 1].name;
                            if (rules[0].right[i].first.size() -1 == get<1>(digram_position[digram_to_erase]))
                                digram_position.erase(digram_to_erase);
                        }
                        rules[0].right[i].first.pop_back();
                        rules[0].right[i].first.push_back(non_terminals[get<0>(digram_position[digram_to_look])]);
                        rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].second.first += 1.0;
                    }
                    if (0==1)
                        cout<< "Found DigramPosition saux2." << endl;
                    /*check_digram_position_integrity_by_rules(digram_position);
                    check_digram_position_integrity(digram_position);*/
                    enforce_rule_utility(digram_position, i);
                    /*check_digram_position_integrity_by_rules(digram_position);
                    check_digram_position_integrity(digram_position);*/

                    if (rules[0].right[i].first.size() > 1)
                        digram = rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " +rules[0].right[i].first[rules[0].right[i].first.size() - 1].name;
                    else
                       digram = "";

                    if (digram_position.find("") != digram_position.end()) {
                        cout << "Achou";
                    }
                }
                if (rules[0].right[i].first.size() > 1) {
                    digram_map.insert(make_pair(digram, 0));
                    digram_position.insert(make_pair(digram, make_tuple(0,rules[0].right[i].first.size() - 1, i)));
                } else
                    digram = "";
            }
            else
                dicard++;


        }
        //enforce_rule_utility();
        //check_digram_position_integrity(digram_position);
        //check_digram_position_integrity_by_rules(digram_position);
    }

}

bool Grammar::Grammar::enforce_digram_uniqueness(std::unordered_map<std::string, int> & digram_map, std::unordered_map<std::string, std::tuple<int, int, int>> & digram_position, int i) {
    digram_map[rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " +rules[0].right[i].first[rules[0].right[i].first.size() - 1].name] = n_non_terminals;
    n_non_terminals++;
    Symbol::Symbol nnt = Symbol::Symbol("NT" + to_string(n_non_terminals-1), n_non_terminals - 1, false, false);
    non_terminals.push_back(nnt);
    vector<Symbol::Symbol> lhs;
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;
    lhs.push_back(nnt);
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> rhss;
    rhs.first.push_back(Symbol::Symbol(rules[0].right[i].first[rules[0].right[i].first.size() - 2]));
    rhs.first.push_back(Symbol::Symbol(rules[0].right[i].first[rules[0].right[i].first.size() - 1]));
    rhs.second.first = 2.0;
    rhs.second.second = ALFA;
    rhss.push_back(rhs);
    rules.push_back(Rule::Rule(lhs, rhss));
    string digram_to_look = rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " +rules[0].right[i].first[rules[0].right[i].first.size() - 1].name;
    string backward_digram;
    string backward_digram_erase;


    if (get<1>(digram_position[digram_to_look]) >1) {
        backward_digram = rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first[get<1>(digram_position[digram_to_look])-2].name + " " + lhs[0].name;
        backward_digram_erase = rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first[get<1>(digram_position[digram_to_look])-2].name + " " + rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first[get<1>(digram_position[digram_to_look])-1].name;
    }
    string forward_digram;
    string forward_digram_erase;
    if (get<1>(digram_position[digram_to_look])+1 < rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.size()) {
        forward_digram = lhs[0].name + " " + rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first[get<1>(digram_position[digram_to_look])+1].name;
        forward_digram_erase  = rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first[get<1>(digram_position[digram_to_look])].name + " " + rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first[get<1>(digram_position[digram_to_look])+1].name;
    }
    rules[0].right[i].first.pop_back();

    // inicio inserção nova regra no lugar do digrama repetido
    if (get<1>(digram_position[digram_to_look]) == 0)
        rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.insert(rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.begin()  +get<1>(digram_position[digram_to_look]), lhs.begin(), lhs.end());
    else
        rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.insert(rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.begin()  +get<1>(digram_position[digram_to_look])-1, lhs.begin(), lhs.end());
    rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.erase(rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.begin()  +get<1>(digram_position[digram_to_look]));
    rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.erase(rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.begin()  +get<1>(digram_position[digram_to_look]));
    // fim inserção nova regra no lugar do digrama repetido

    //digram_map.erase(rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " + rules[0].right[0].first[rules[0].right[0].first.size() - 1].name);

    if (rules[0].right[i].first.size() > 1) {
        string saux = rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " + rules[0].right[i].first[rules[0].right[i].first.size() - 1].name;
        if (rules[0].right[i].first[rules[0].right[i].first.size() - 2].name.compare(rules[0].right[i].first[rules[0].right[i].first.size() - 1].name) == 0) {
            if( rules[0].right[i].first.size() >2) {
                if (rules[0].right[i].first[rules[0].right[i].first.size() - 3].name.compare(rules[0].right[i].first[rules[0].right[i].first.size() - 2].name) == 0) {
                    cout << "digramas consecutivos iguais" << endl;
                } else
                    digram_position.erase(rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " + rules[0].right[i].first[rules[0].right[i].first.size() - 1].name);
            }
            else
                digram_position.erase(rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " + rules[0].right[i].first[rules[0].right[i].first.size() - 1].name);
        }
        else
            digram_position.erase(rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " + rules[0].right[i].first[rules[0].right[i].first.size() - 1].name);

    }
    if (get<1>(digram_position[digram_to_look])+1 < rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.size()) {
        if (digram_position.find(forward_digram_erase) != digram_position.end()) {
            for (auto m: digram_position) {
                if (get<0>(m.second) == get<0>(digram_position[forward_digram_erase]))
                    if (get<2>(m.second) == get<2>(digram_position[forward_digram_erase]))
                        if (get<1>(m.second) > get<1>(digram_position[forward_digram_erase])) {
                            get<1>(digram_position[m.first])--;
                        }
            }
        }
    }
    rules[0].right[i].first.pop_back();
    rules[0].right[i].first.push_back(nnt);
    if (get<1>(digram_position[digram_to_look]) >1) {
        digram_map.insert(make_pair(backward_digram, 0));
        digram_position.insert(make_pair(backward_digram, make_tuple(get<0>(digram_position[digram_to_look]), get<1>(digram_position[digram_to_look])-1, get<2>(digram_position[digram_to_look]))));
        if (get<1>(digram_position[digram_to_look])-1 == get<1>(digram_position[backward_digram_erase])) {
            digram_map.erase(backward_digram_erase);
            digram_position.erase(backward_digram_erase);
        }

    }
    forward_digram = lhs[0].name + " " + rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first[get<1>(digram_position[digram_to_look])].name;
    //TODO atenção com esse +1 (ou falta dele) antes do <, talvez tenha bug.
    if (get<1>(digram_position[digram_to_look]) < rules[get<0>(digram_position[digram_to_look])].right[get<2>(digram_position[digram_to_look])].first.size()){
        digram_map.insert(make_pair(forward_digram, get<2>(digram_position[digram_to_look])));
        digram_position.insert(make_pair(forward_digram, digram_position[digram_to_look]));
    }

    if (digram_position.find(forward_digram_erase) != digram_position.end()) {
        string possibleTripleDigram = rules[get<0>(digram_position[forward_digram_erase])].
                                  right[get<2>(digram_position[forward_digram_erase])].
                                  first[get<1>(digram_position[forward_digram_erase])-1].name + " " +
                                  rules[get<0>(digram_position[forward_digram_erase])].
                                  right[get<2>(digram_position[forward_digram_erase])].
                                  first[get<1>(digram_position[forward_digram_erase])].name;
        if (forward_digram_erase.compare(possibleTripleDigram) != 0)
            digram_position.erase(forward_digram_erase);
    }

    //digram_position.erase(digram_to_look);
    get<0>(digram_position[digram_to_look])= nnt.id;
    get<2>(digram_position[digram_to_look]) = 0;
    get<1>(digram_position[digram_to_look]) = 1;
    return true;
}



bool Grammar::Grammar::verify_duplicate_digram(std::unordered_map<std::string, std::tuple<int, int, int>> digram_map, string digram) {
    if (digram_map.find(digram) == digram_map.end()) {
        return false;
    }
    else return true;
}



/*
void Grammar::Grammar::prob_sequitur() {
    non_terminals.clear();
    rules.clear();
    n_non_terminals = 0;
    n_non_terminals = 1;
    Symbol::Symbol start = Symbol::Symbol("NT"+ to_string(n_non_terminals-1), n_non_terminals-1, false, false);
    non_terminals.push_back(start);
    vector<Symbol::Symbol> lhs;
    lhs.push_back(start);
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> rhss;
    rules.push_back(Rule::Rule(lhs,rhss));
    for (int i = 0; i < words.size(); i++) {
        auto w = words[i];
        std::unordered_map<string, int> digram_map;
        //td::unordered_map<string, pair<int, int>> digram_position;
        std::unordered_map<string, tuple<int, int, int>> digram_position;
        std::unordered_map<int, int> rule_map;
        std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;
        rhs.first.push_back(Symbol::Symbol(w[0]));
        rhs.second.first = 1.0;
        rhs.second.second = ALFA;
        rhss.push_back(rhs);
        rules[0].right.push_back(rhs);
        int dicard = 0;
        for (auto s: w) {
            print_grammar();
            if (!(dicard <1)) {
                rules[0].right[0].first.push_back(s);
                string digram = rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " +rules[0].right[0].first[rules[0].right[0].first.size() - 1].name;
                while (verify_duplicate_digram(digram_position, digram)) {
                    if (!((get<0>(digram_position[digram]) != 0) && (rules[get<0>(digram_position[digram])].right[0].first.size() == 2 ))) {
                        if (rules[0].right[0].first.size() - 1 - get<1>(digram_position[digram]) > 1) {
                            if (!rules[0].right[0].first[rules[0].right[0].first.size() - 2].terminal)
                                rules[rules[0].right[0].first[rules[0].right[0].first.size() - 2].id].right[0].second.first -= 1.0;
                            if (!rules[0].right[0].first[rules[0].right[0].first.size() - 1].terminal)
                                rules[rules[0].right[0].first[rules[0].right[0].first.size() - 1].id].right[0].second.first -= 1.0;
                            enforce_digram_uniqueness(digram_map,  digram_position);
                        }
                        else
                            break;
                    } else {
                        string digram_to_look = rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " + rules[0].right[0].first[rules[0].right[0].first.size() - 1].name;
                        if (!rules[0].right[0].first[rules[0].right[0].first.size() - 2].terminal)
                            rules[rules[0].right[0].first[rules[0].right[0].first.size() - 2].id].right[0].second.first -= 1.0;
                        if (!rules[0].right[0].first[rules[0].right[0].first.size() - 1].terminal)
                            rules[rules[0].right[0].first[rules[0].right[0].first.size() - 1].id].right[0].second.first -= 1.0;
                        rules[0].right[0].first.pop_back();
                        digram_map.erase(rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " + rules[0].right[0].first[rules[0].right[0].first.size() - 1].name);digram_position.erase(rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " + rules[0].right[0].first[rules[0].right[0].first.size() - 1].name);
                        rules[0].right[0].first.pop_back();
                        rules[0].right[0].first.push_back(non_terminals[get<0>(digram_position[digram_to_look])]);
                        rules[get<0>(digram_position[digram_to_look])].right[0].second.first += 1.0;
                    }
                    digram = rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " +rules[0].right[0].first[rules[0].right[0].first.size() - 1].name;
                    enforce_rule_utility(digram_position);
                }
                digram_map.insert(make_pair(digram, 0));
                digram_position.insert(make_pair(digram, make_tuple(0,rules[0].right[0].first.size() - 1, i)));
            }
            else
                dicard++;
        }
        //enforce_rule_utility();
    }
}

*/


/*
bool Grammar::Grammar::enforce_digram_uniqueness(std::unordered_map<std::string, int> & digram_map, std::unordered_map<std::string, std::tuple<int, int, int>> & digram_position, int i) {
    digram_map[rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " +rules[0].right[0].first[rules[0].right[0].first.size() - 1].name] = n_non_terminals;
    n_non_terminals++;
    Symbol::Symbol nnt = Symbol::Symbol("NT" + to_string(n_non_terminals-1), n_non_terminals - 1, false, false);
    non_terminals.push_back(nnt);
    vector<Symbol::Symbol> lhs;
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;
    lhs.push_back(nnt);
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> rhss;
    rhs.first.push_back(Symbol::Symbol(rules[0].right[0].first[rules[0].right[0].first.size() - 2]));
    rhs.first.push_back(Symbol::Symbol(rules[0].right[0].first[rules[0].right[0].first.size() - 1]));
    rhs.second.first = 2.0;
    rhs.second.second = ALFA;
    rhss.push_back(rhs);
    rules.push_back(Rule::Rule(lhs, rhss));
    string digram_to_look = rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " +rules[0].right[0].first[rules[0].right[0].first.size() - 1].name;
    string backward_digram;
    string backward_digram_erase;
    if (get<1>(digram_position[digram_to_look]) >1) {
        */
/*backward_digram = rules[0].right[0].first[digram_position[digram_to_look].second-2].name + " " + lhs[0].name;
        backward_digram_erase = rules[0].right[0].first[digram_position[digram_to_look].second-2].name + " " + rules[0].right[0].first[digram_position[digram_to_look].second-1].name;*//*

        backward_digram = rules[get<0>(digram_position[digram_to_look])].right[0].first[get<1>(digram_position[digram_to_look])-2].name + " " + lhs[0].name;
        backward_digram_erase = rules[get<0>(digram_position[digram_to_look])].right[0].first[get<1>(digram_position[digram_to_look])-2].name + " " + rules[get<0>(digram_position[digram_to_look])].right[0].first[get<1>(digram_position[digram_to_look])-1].name;
    }
    //novo forward digram
    string forward_digram = lhs[0].name + " " + rules[get<0>(digram_position[digram_to_look])].right[0].first[get<1>(digram_position[digram_to_look])+1].name;
    string forward_digram_erase  = rules[get<0>(digram_position[digram_to_look])].right[0].first[get<1>(digram_position[digram_to_look])].name + " " + rules[get<0>(digram_position[digram_to_look])].right[0].first[get<1>(digram_position[digram_to_look])+1].name;
    rules[0].right[0].first.pop_back();
    rules[get<0>(digram_position[digram_to_look])].right[0].first.insert(rules[get<0>(digram_position[digram_to_look])].right[0].first.begin()  +get<1>(digram_position[digram_to_look])-1, lhs.begin(), lhs.end());
    rules[get<0>(digram_position[digram_to_look])].right[0].first.erase(rules[get<0>(digram_position[digram_to_look])].right[0].first.begin()  +get<1>(digram_position[digram_to_look]));
    rules[get<0>(digram_position[digram_to_look])].right[0].first.erase(rules[get<0>(digram_position[digram_to_look])].right[0].first.begin()  +get<1>(digram_position[digram_to_look]));
    digram_map.erase(rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " + rules[0].right[0].first[rules[0].right[0].first.size() - 1].name);
    digram_position.erase(rules[0].right[0].first[rules[0].right[0].first.size() - 2].name + " " + rules[0].right[0].first[rules[0].right[0].first.size() - 1].name);
    if (get<1>(digram_position[digram_to_look]) < rules[get<0>(digram_position[digram_to_look])].right[0].first.size()) {
        for (auto m: digram_position) {
            if (get<0>(m.second) == get<0>(digram_position[forward_digram_erase]))
                if (get<1>(m.second) > get<1>(digram_position[forward_digram_erase])) {
                    get<1>(digram_position[m.first])--;
                }

        }
    }
    rules[0].right[0].first.pop_back();
    rules[0].right[0].first.push_back(nnt);
    if (get<1>(digram_position[digram_to_look]) >1) {
        digram_map.erase(backward_digram_erase);
        digram_position.erase(backward_digram_erase);
        digram_map.insert(make_pair(backward_digram, 0));
        digram_position.insert(make_pair(backward_digram, make_tuple(0, get<1>(digram_position[digram_to_look])-1, 0)));
    }
    forward_digram = lhs[0].name + " " + rules[0].right[0].first[get<1>(digram_position[digram_to_look])].name;
    forward_digram = lhs[0].name + " " + rules[get<0>(digram_position[digram_to_look])].right[0].first[get<1>(digram_position[digram_to_look])].name;
    if (get<1>(digram_position[digram_to_look]) < rules[get<0>(digram_position[digram_to_look])].right[0].first.size()){
        digram_map.insert(make_pair(forward_digram, 0));
        digram_position.insert(make_pair(forward_digram, digram_position[digram_to_look]));
    }
    digram_map.erase(forward_digram_erase);
    digram_position.erase(forward_digram_erase);
    //digram_position.erase(digram_to_look);
    get<0>(digram_position[digram_to_look])= nnt.id;
    get<1>(digram_position[digram_to_look]) = 0;
    return true;
}
*/


void Grammar::Grammar::enforce_rule_utility(std::unordered_map<std::string, std::tuple<int, int, int>> & digram_position, int ii) {
    std::vector<Rule::Rule>::iterator itR;
    int cont = 0;
    for (itR = rules.begin();  itR != rules.end(); itR++) {
        Rule::Rule r = (*itR);
        if (r.right[0].second.first == 1.0 && r.left[0].id != 0) {
            cont++;
            if (r.left[0].id == 288)
                cout << "found" << endl;
            if(r.left[0].id == 289)
                cout << "found" << endl;
            if(r.left[0].id == 287)
                cout << "found" << endl;
            string saux = "NT174 NT174";
            (*itR).right[0].second.first = 0.0;
            std::vector<Rule::Rule>::iterator itRule;
            for (itRule = rules.begin();  itRule != rules.end(); itRule++) {
                int index = 0;
                if (itRule == rules.begin())
                    index = ii;
                bool rule_found = false;
                for (int i = 0; i < (*itRule).right[index].first.size(); i++) {
                    if (!(*itRule).right[index].first[i].terminal) {
                        if (r.left[0].id == (*itRule).right[index].first[i].id) {

                            string digram = r.right[0].first[0].name + " "+ r.right[0].first[1].name;
                            if( (*itRule).left[0].id == 288)
                                cout << "Achou" << endl;
                            //string digram = r.right[index].first[0].name + " "+ r.right[index].first[1].name;
                            //  troquei inder por 0 na regra r.right[]
                            (*itRule).right[index].first.insert((*itRule).right[index].first.begin()+i, r.right[0].first.begin(), r.right[0].first.end());
                            (*itRule).right[index].first.erase((*itRule).right[index].first.begin()+i+r.right[0].first.size());
                            string previous_digram = "";
                            for (int j = 1; j < (*itRule).right[0].first.size(); j++) {
                                string new_digram_aux = (*itRule).right[0].first[j-1].name + " " + (*itRule).right[0].first[j].name;
                                if (digram_position.find(new_digram_aux) != digram_position.end()) {
                                    if (get<2>(digram_position[new_digram_aux]) != 0 ) {
                                        string t1 = "G:maj C:min";
                                        t1 = "C:min G:min";
                                        t1 = "G:min C:min";
                                        //print_grammar();
                                        cout << "Achou" << endl;
                                    }
                                    if (new_digram_aux.compare(previous_digram) != 0) {
                                        get<0>(digram_position[new_digram_aux]) = (*itRule).left[0].id;
                                        get<1>(digram_position[new_digram_aux]) = j;
                                        get<2>(digram_position[new_digram_aux]) = index;
                                    }
                                    previous_digram = new_digram_aux;
                                }
                                else
                                    digram_position.insert(make_pair(new_digram_aux, make_tuple( (*itRule).left[0].id,  j, 0)));

                            }
                            get<0>(digram_position[digram]) = (*itRule).left[0].id;
                            get<1>(digram_position[digram]) = i+1;
                            //troquei 0 por index nos inserts abaixo
                            if (i > 0) {
                                digram = (*itRule).right[index].first[i-1].name + " " + r.right[0].first[0].name;
                                digram_position.insert(make_pair(digram, make_tuple( (*itRule).left[0].id, i, index)));
                                digram_position.erase((*itRule).right[index].first[i-1].name+ " " + r.left[0].name );
                            } else if (i < (*itRule).right[index].first.size()-1 ){
                                digram = r.right[0].first[r.right[0].first.size()-1].name + " " + (*itRule).right[index].first[i + r.right[0].first.size()].name;
                                digram_position.insert(make_pair(digram, make_tuple( (*itRule).left[0].id,  i+ r.right[0].first.size(), index)));
                                digram_position.erase(r.left[0].name+ " " + (*itRule).right[index].first[i + r.right[0].first.size()].name );
                            }
                            rule_found = true;
                            break;
                        }
                    }
                }
                if(rule_found)
                    break;
            }

        }
    }
    //cout << "regras removidas: " << cont << endl;
}
void Grammar::Grammar::convert_to_cnf() {
    //IMPORTANTE: ao testar se a gramática na cnf gera a base de dados, comentar o trecho abaixo qm que regras iniciais iguais são agrupadas.
    //GROUP EQUAL INITIAl RULES
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    for (int j = 0; j < rules[0].right.size(); j++) {
    //for (itRight = rules[0].right.begin(); itRight != rules[0].right.end(); itRight++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight2;
        //for (itRight2 = rules[0].right.begin()+1; itRight2 != rules[0].right.end(); itRight2++) {
        for (int i = j+1; i < rules[0].right.size(); i++) {
            if (equal_rhs(rules[0].right[j].first,  rules[0].right[i].first)) {
                rules[0].right[j].second.first += 1.0;
                rules[0].right.erase(rules[0].right.begin()+i);
                i--;
            }
        }
    }

    //TERMINAL RULES
    vector<Rule::Rule> rulesTerm;
    for (auto t: terminals) {
        n_non_terminals++;
        Symbol::Symbol nnt = Symbol::Symbol("NT" + to_string(n_non_terminals-1), n_non_terminals - 1, false, false);
        non_terminals.push_back(nnt);
        vector<Symbol::Symbol> lTerminal;
        lTerminal.push_back(nnt);
        vector<Symbol::Symbol> rTerminal;
        rTerminal.push_back(t);
        vector<pair<vector<Symbol::Symbol>,pair<double, double>>> rhsTerminal;
        rhsTerminal.push_back(make_pair(rTerminal, make_pair(0.0,0.1)));
        Rule::Rule r = Rule::Rule(lTerminal, rhsTerminal);
        rulesTerm.push_back(r);
    }
    //print_grammar();
    //TERM
    vector<Rule::Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            for (int i = 0; i < (*itRight).first.size(); i++) {
                if ((*itRight).first[i].terminal) {
                    vector<pair<std::vector<Symbol::Symbol>, pair<double, double>>> rTerminal;
                    (*itRight).first[i] = Symbol::Symbol(rulesTerm[(*itRight).first[i].id].left[0]);
                    //mudei linha abaixo
                    rulesTerm[(*itRight).first[i].id - rulesTerm[0].left[0].id].right[0].second.first += (*itRight).second.first;
                }
            }
        }
    }
    rules.insert(rules.end(), rulesTerm.begin(), rulesTerm.end());
    //print_grammar();
    rulesTerm.clear();

    //BIN (testar com linguagens com palavras maiores que 2)
    vector<Rule::Rule> rulesAux;
    rulesAux.insert(rulesAux.begin(), rules.begin(), rules.end());
    for (itRule = rulesAux.begin(); itRule != rulesAux.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            while ((*itRight).first.size() > 2) {
                vector<Symbol::Symbol> rhs_cnf((*itRight).first.begin()+1, (*itRight).first.begin()+3);
                vector<Symbol::Symbol> lTerminal;
                vector<pair<vector<Symbol::Symbol>,pair<double, double>>> rhsTerminal;
                int nt_index = verify_rule_existence_cnf(rhs_cnf, rulesTerm);
                Rule::Rule r = Rule::Rule(lTerminal,rhsTerminal);
                if (nt_index != -1) {
                    r = rulesTerm[nt_index - rules.size()];
                    rulesTerm[nt_index - rules.size()].right[0].second.first += 1.0;
                } else {
                    vector<Symbol::Symbol> lTerminal;
                    n_non_terminals++;
                    Symbol::Symbol nnt = Symbol::Symbol("NT" + to_string(n_non_terminals-1), n_non_terminals - 1, false, false);
                    non_terminals.push_back(nnt);
                    lTerminal.push_back(nnt);
                    rhsTerminal.push_back(make_pair(rhs_cnf, make_pair(1.0,0.1)));
                    r = Rule::Rule(lTerminal, rhsTerminal);
                    rulesTerm.push_back(r);
                }
                (*itRight).first[1] = Symbol::Symbol(r.left[0]);
                (*itRight).first.erase((*itRight).first.begin()+2);
            }
        }
    }
    rules = rulesAux;
    rules.insert(rules.end(), rulesTerm.begin(), rulesTerm.end());
    rulesTerm.clear();

    //UNIT
    for (itRule = rules.begin(); itRule != rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>::iterator itRight;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            if ((*itRight).first.size() == 1) {
                if (!(*itRight).first[0].terminal) {
                    int iRule = (*itRight).first[0].id;
                    (*itRight).first.clear();
                    (*itRight).first.insert((*itRight).first.begin(), rules[iRule].right[0].first.begin(), rules[iRule].right[0].first.end());
                    rules[iRule].right[0].second.first -= (*itRight).second.first;
                }
            }
        }
    }
    //print_grammar();
}
int Grammar::Grammar::verify_rule_existence_cnf(std::vector<Symbol::Symbol> rhs_cnf, std::vector<Rule::Rule> rules_nt_cnf) {
    int flag = -1;
    for (auto r: rules_nt_cnf) {
        if (rhs_cnf.size() != r.right[0].first.size())
            return -1;
        flag = r.left[0].id;
        for (int i = 0; i < rhs_cnf.size(); i++) {
            if (!rhs_cnf[i].equal_symbol(r.right[0].first[i]))
                return -1;
        }
    }
    return flag;
}

bool Grammar::Grammar::equal_rhs(std::vector<Symbol::Symbol> rh1, std::vector<Symbol::Symbol> rh2) {
    if (rh1.size() != rh2.size())
        return false;
    for (int i = 0; i < rh1.size(); i++) {
        if (!rh1[i].equal_symbol(rh2[i]))
            return false;
    }
    return true;
}
//TODO fazer algoritmo de contagem de bombeamentos por slices igual no caderno
//TODO fazer programação dinâmica constexpr
void Grammar::Grammar::count_pumping_str(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, int sub_amount, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word) {

    tuple<vector<Symbol::Symbol>, vector<Symbol::Symbol>, vector<Symbol::Symbol>, vector<Symbol::Symbol>, vector<Symbol::Symbol>> pumping_str;

    switch (word.size()) {
        case 1:
            count_pumping_size1(map, word, map_pump_to_word, original_word); break;
        case 2:
            count_pumping_size2(map, word, map_pump_to_word, original_word); break;
        case 3:
            count_pumping_size3(map, word, map_pump_to_word, original_word); break;
        case 4:
            count_pumping_size4(map, word, map_pump_to_word, original_word); break;
        case 5:
            count_pumping_size5(map, word, map_pump_to_word, original_word); break;
    }
    for (int i = sub_amount; i <word.size()-1 ; i ++) {
        vector<Symbol::Symbol> aux_word;
        aux_word.insert(aux_word.end(), word.begin(), word.end());
        aux_word[i].name += " "+aux_word[i+1].name;
        aux_word.erase(aux_word.begin()+i+1);
        count_pumping_str(map, aux_word, i, map_pump_to_word, original_word);
    }
}

void Grammar::Grammar::count_pumping_size1(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word) {
    Symbol::Symbol word_0_pumped = verify_pumping_times(word[0]);
    string pumping;
    if (word_0_pumped.name.compare("")  != 0){
        pumping = "|" +word_0_pumped.name + "|||"; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
        /*pumping = "|||" +word_0_pumped.name + "|"; map[pumping] +=1;*/
    }
}
/*
void Grammar::Grammar::count_pumping_size2(unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word) {
    string pumping = "|" +word[1].name + "|||"; map[pumping] +=1;
    pumping =  "|||" +word[1].name + "|"; map[pumping] +=1;
    pumping = "|" + word[0].name + "|||" ; map[pumping] +=1;
    pumping = "|" + word[0].name + "||" +word[1].name + "|"; map[pumping] +=1;
    pumping = "|" + word[0].name + "|||"; map[pumping] +=1;
    pumping = "|||"  + word[1].name + "|"; map[pumping] +=1;
    pumping = "|||" + word[0].name + "|" ; map[pumping] +=1;

}

void Grammar::Grammar::count_pumping_size3(unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word) {
    string pumping =  "|" +word[1].name + "|" + "||"; map[pumping] +=1;
    pumping = "|" +word[1].name +"||" +word[2].name +"|"; map[pumping] +=1;
    pumping = "|" +word[1].name +"|||" ; map[pumping] +=1;
    pumping = "|" + word[0].name + "||" +word[2].name + "|"; map[pumping] +=1;
    pumping = "|" + word[0].name + "||" +word[1].name + + "|"; map[pumping] +=1;
    pumping = "|" + word[0].name + "|||"; map[pumping] +=1;
    pumping = "|||" + word[1].name + "|"; map[pumping] +=1;
}

void Grammar::Grammar::count_pumping_size4(unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word) {
    string pumping =  "|" +word[1].name + "||" +word[3].name +"|"; map[pumping] +=1;
    pumping =  "|" +word[1].name + "|||"; map[pumping] +=1;
    pumping =  "|" +word[1].name + "||" +word[2].name + "|"; map[pumping] +=1;
    pumping = "||" +word[1].name + "|" +word[2].name + "|" ; map[pumping] +=1;
    pumping = "|" + word[0].name + "||" +word[2].name + "|" ; map[pumping] +=1;
}
*/


void Grammar::Grammar::count_pumping_size2(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word) {
    Symbol::Symbol word_0_pumped = verify_pumping_times(word[0]);
    Symbol::Symbol  word_1_pumped = verify_pumping_times(word[1]);
    string pumping;
    if (word_0_pumped.name.compare("")  != 0){
        pumping = "|" + word_0_pumped.name + "|||" +word[1].name ; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
        /*pumping = "|" + word_0_pumped.name + "|" + word[1].name + "||"; map[pumping] +=1;
        pumping = "|||" + word_0_pumped.name + "|" + word[1].name ; map[pumping] +=1;*/
    } else if (word_1_pumped.name.compare("")  != 0) {
        pumping = word[0].name + "|" +word_1_pumped.name + "|||"; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
        /*pumping = word[0].name + "|||" +word_1_pumped.name + "|"; map[pumping] +=1;
        pumping = "||" + word[0].name + "|" + word_1_pumped.name + "|"; map[pumping] +=1;*/
    } else if (word_0_pumped.name.compare("")  != 0 && word_1_pumped.name.compare("")) {
        pumping = "|" + word_0_pumped.name + "||" +word_1_pumped.name + "|"; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    }
}

void Grammar::Grammar::count_pumping_size3(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word) {
    string pumping;
    Symbol::Symbol word_0_pumped = verify_pumping_times(word[0]);
    Symbol::Symbol word_1_pumped = verify_pumping_times(word[1]);
    Symbol::Symbol word_2_pumped = verify_pumping_times(word[2]);
    if (word_1_pumped.name.compare("")  != 0) {
        pumping = word[0].name + "|" +word_1_pumped.name + "|" +word[2].name + "||"; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
        /*pumping = word[0].name + "|" +word_1_pumped.name +"|||" +word[2].name; map[pumping] +=1;
        pumping = "||" + word[0].name + "|" + word_1_pumped.name + "|" +word[2].name ; map[pumping] +=1;
        pumping = word[0].name + "|||"  + word_1_pumped.name + "|" +word[2].name ; map[pumping] +=1;*/
    } else if (word_1_pumped.name.compare("")  != 0 && word_2_pumped.name.compare("")  != 0) {
        pumping = word[0].name + "|" +word_1_pumped.name +"||" +word_2_pumped.name +"|"; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    } else if (word_0_pumped.name.compare("") != 0 && word_1_pumped.name.compare("")  != 0) {
        pumping = "|" + word_0_pumped.name + "||" + word_1_pumped.name + +"|" + word[2].name; map[pumping] += 1;
        map_pump_to_word[pumping].push_back(original_word);
    } else if(word_0_pumped.name.compare("") != 0) {
        pumping = "|" + word_0_pumped.name + "|" + word[1].name + "||"  +word[2].name; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    } else if (word_0_pumped.name.compare("") != 0 && word_2_pumped.name.compare("")  != 0) {
        pumping = "|" + word_0_pumped.name + "|" +word[1].name + + "|" +word_2_pumped.name + "|"; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    }
}

void Grammar::Grammar::count_pumping_size4(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word) {
    string pumping;
    Symbol::Symbol word_0_pumped = verify_pumping_times(word[0]);
    Symbol::Symbol word_1_pumped = verify_pumping_times(word[1]);
    Symbol::Symbol word_2_pumped = verify_pumping_times(word[2]);
    Symbol::Symbol word_3_pumped = verify_pumping_times(word[3]);
    if (word_1_pumped.name.compare("")  != 0) {
        pumping = word[0].name + "|" +word_1_pumped.name + "|" +word[2].name + "||" +word[3].name ; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    } else if (word_1_pumped.name.compare("")  != 0 && word_3_pumped.name.compare("")  != 0) {
        pumping = word[0].name + "|" +word[1].name + "|" +word[2].name + "|" +word_3_pumped.name +"|"; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    } else if (word_1_pumped.name.compare("")  != 0 && word_2_pumped.name.compare("")  != 0) {
        pumping = word[0].name + "|" +word_1_pumped.name + "||" +word_2_pumped.name + "|" +word[3].name; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    } else if (word_2_pumped.name.compare("")  != 0) {
        pumping = word[0].name + "||" +word[1].name + "|" +word_2_pumped.name + "|" +word[3].name; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    } else if (word_0_pumped.name.compare("")  != 0 && word_2_pumped.name.compare("")  != 0) {
        pumping = "|" + word_0_pumped.name + "|" +word[1].name + "|" +word_2_pumped.name + "|" +word[3].name; map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    }
}
void Grammar::Grammar::count_pumping_size5(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word, std::vector<Symbol::Symbol> &original_word) {
    Symbol::Symbol word_1_pumped = verify_pumping_times(word[1]);
    Symbol::Symbol word_3_pumped = verify_pumping_times(word[3]);
    if (word_1_pumped.name.compare("")  != 0 && word_3_pumped.name.compare("")  != 0) {
        string pumping = word[0].name +"|" +word_1_pumped.name + "|"+word[2].name +"|" + word_3_pumped.name +"|"+word[4].name;
        map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    }

}

Symbol::Symbol Grammar::Grammar::verify_pumping_times(Symbol::Symbol s) {
    stringstream check1(s.name);
    string intermediate;
    vector<string> tokens;
    while(getline(check1, intermediate, ' '))
            tokens.push_back(intermediate);
    bool flag = true;
    if(tokens.size()<2)
        flag = false;
    for (int i = 1; i < tokens.size()/2+1; i++) {
        if (tokens.size()%i == 0) {
            flag = true;
            int j = 0;
            for (j = 0; j < tokens.size()-i; j+=i) {
                for (int k = j; k <j+i; k++) {
                    if (tokens[k].compare(tokens[k+i]) != 0) {
                        flag = false;
                        break;
                    }
                }
                if(!flag)
                    break;

            }

            if (flag == true) {
                if (j < i ) {
                    s.name = "";
                    return s;
                }
                s.name = "";
                for (int k = 0; k < i; k++)
                    s.name += tokens[k] + " ";
                s.name = s.name.substr(0, s.name.size()-1);
                return s;
            }

        }
    }
    if (flag == false)
        s.name = "";
    return s;
}

bool sortbysec(const pair<string ,int> &a,
               const pair<string,int> &b)
{
    return (a.second > b.second);
}

void Grammar::Grammar::pumping_inference(unordered_map<std::string, int> &map, std::unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> & map_pump_to_word) {
    //TODO: Reduzir o count na regra do NT
    vector<bool> pumped_rule;
    for (int i = 0; i < rules[0].right.size(); i++)
        pumped_rule.push_back(false);
    vector<pair<string, int>> ordered;
    for (auto w : map)
        ordered.push_back(w);
    std::sort(ordered.begin(), ordered.end(), sortbysec);
    for (auto s: ordered) {
        if (s.second <5)
            break;
        cout << s.first << ": " << s.second << endl;
    }


    unordered_map<string, Symbol::Symbol> non_terminal_string_map;
    int repeated = 0;
    for (int i = 1; i < rules.size(); i++) {
        vector<Symbol::Symbol> yielded = yield_string(rules[i].left);
        string key = "";
        for (int k = 0; k < yielded.size(); k++)
            key += yielded[k].name + " ";
        key = key.substr(0, key.size()-1);
        if (non_terminal_string_map.find(key) == non_terminal_string_map.end())
            non_terminal_string_map[key] = rules[i].left[0];
        else
            cout << "rule: "<< i << " repeated key: " << key << endl;
    }
    vector<pair<vector<Symbol::Symbol>,pair<double, double>>> new_start_right;
    int countRules = 0;
    for (auto u: ordered) {
        //criterio de parada com count de bombeamento menor que 5
        if (u.second < 5)
            break;
        vector<pair<vector<Symbol::Symbol>,pair<double, double>>>::iterator itRight;
        int v;
       /* if (n_non_terminals ==28)
            cout << "testar Aqui. nnt = : " << n_non_terminals << endl;*/
        find_pumping_rules(v, u.first, non_terminal_string_map);
        for (itRight = rules[0].right.begin(); itRight != rules[0].right.end(); itRight++) {
            if (pumped_rule[itRight-rules[0].right.begin()] == false) {
                if (compare_pumping_use(itRight->first, map_pump_to_word[u.first])) {
                    countRules++;
                    //print_grammar();
                    cout << "Pumpumg Rules used: " << countRules << endl;
                    pair<vector<Symbol::Symbol>,pair<double, double>> new_rhs;
                    calculate_new_rule_from_starting_symbol(new_rhs, v, u.first, non_terminal_string_map);
                    int existRhsFlag = false;
                    for (auto & r: new_start_right) {
                        if (equal_word(r.first, new_rhs.first)) {
                            pumped_rule[itRight-rules[0].right.begin()] = true;
                            r.second.first += itRight->second.first;
                            rules[v].right[0].second.first += itRight->second.first;
                            //rules[v].right[1].second.first += itRight->second.first;
                            itRight->second.first = 0;
                            existRhsFlag = true;
                            break;
                        }
                    }
                    if (!existRhsFlag) {
                        new_rhs.second.first = itRight->second.first;
                        //rules[v].right[0].second.first += itRight->second.first;
                        //TO DO: ajustar a linha abaixo pra quantidade de vezes que é gerada a subpalavra bombeada
                        //rules[v].right[0].second.first += itRight->second.first;
                        rules[v].right[0].second.first += 0.5;
                        //Se pertmir bombeamentos com apenas 1 utilização, vai dar merda
                        rules[v].right[1].second.first += itRight->second.first-0.5;
                        itRight->second.first = 0;
                        pumped_rule[itRight-rules[0].right.begin()] = true;
                        new_start_right.push_back(new_rhs);
                    }
                    for (auto r: itRight->first)
                        if (r.id != v)
                            rules[r.id].right[0].second.first -= new_rhs.second.first;
                    /*for (auto r: new_rhs.first) {
                        if (r.id != v)
                            rules[r.id].right[0].second.first += new_rhs.second.first;
                    }*/




                }
            }
        }
    }
    rules[0].right.insert(rules[0].right.end(), new_start_right.begin(), new_start_right.end());
    print_grammar();

    //para cada bombeamento com contador > que criterio (cobertura de palavras, mínimo de variáveis, %de bombeadas em relação ao tamanho do conjunto de treinamento)
        //para cada regra em NT0
            //criar nova regra do bombeamento
            //derivar palavra e comparar com o uso do bombeamento do mapa
            //se usar bombeamento
                //incrementar contador na nova regra
                //decrementar uso da regra que não será mais usada
                //marcar regra se ela usar o bombeamento.


}
bool Grammar::Grammar::compare_pumping_use(std::vector<Symbol::Symbol> vector, std::vector<std::vector<Symbol::Symbol>> pumping_string) {
    std::vector<Symbol::Symbol> word = yield_string(vector);
    for (auto v: pumping_string) {
        if (equal_word(word, v))
            return true;
    }
    return false;
}
void Grammar::Grammar::find_pumping_rules(int &v, std::string uvxyz, std::unordered_map<std::string, Symbol::Symbol> & non_terminal_string_map) {
    v = -1;
    std::vector<Symbol::Symbol> rhs_to_look;
    size_t init_pipe = uvxyz.find("|", 0);
    size_t next_pipe = uvxyz.find("|", init_pipe+1);
    size_t next_bar = uvxyz.substr(init_pipe, next_pipe-init_pipe).find(" ", 0);
    string pumping_v = "";
    size_t bar = 1;//uvxyz.substr(0, next_pipe-init_pipe).find(" ", 0);
    if(next_pipe-init_pipe<=1)
        bar = string::npos;
    while (bar != string::npos) {
        if (next_bar != string::npos) {
            pumping_v += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, next_bar);
            bar = next_bar;
            next_bar = uvxyz.substr(init_pipe, next_pipe-init_pipe).find(" ", bar+1);
        }
        else {
            if (uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, string::npos)[0] != ' ')
                pumping_v += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, string::npos);
            else
                pumping_v += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar+1, string::npos);
            bar = next_bar;
        }
    }



    init_pipe = next_pipe;
    next_pipe = uvxyz.find("|", init_pipe+1);
    string npumping_x = "";
    next_bar = uvxyz.substr(init_pipe, next_pipe-init_pipe).find(" ", 0);
    bar = 1;
    if(next_pipe-init_pipe<=1)
        bar = string::npos;
    while (bar != string::npos) {
        if (next_bar != string::npos) {
            npumping_x += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, next_bar);
            bar = next_bar;
            next_bar = uvxyz.substr(init_pipe, next_pipe-init_pipe).find(" ", bar+1);
        }
        else {
            npumping_x += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, string::npos);
            bar = next_bar;
        }
    }
    init_pipe = next_pipe;
    next_pipe = uvxyz.find("|", init_pipe+1);
    next_bar = uvxyz.substr(init_pipe, next_pipe-init_pipe).find(" ", 0);
    string pumping_y = "";
    bar = 1;
    if(next_pipe-init_pipe<=1)
        bar = string::npos;
    while (bar != string::npos) {
        if (next_bar != string::npos) {
            pumping_y += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, next_bar);
            bar = next_bar;
            next_bar = uvxyz.substr(init_pipe, next_pipe-init_pipe).find(" ", bar+1);
        }
        else {
            pumping_y += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, string::npos);
            bar = next_bar;
        }
    }

    for (int i = 0; i < rules.size(); i++) {
        if (npumping_x.empty()) {
            if (!pumping_v.empty()) {
                if (!pumping_y.empty()) {
                    if (rules[i].right[0].first[0].equal_symbol(non_terminal_string_map[pumping_v]))
                        if (rules[i].right[0].first[1].equal_symbol(rules[i].left[0]))
                            if (rules[i].right[0].first[2].equal_symbol(non_terminal_string_map[pumping_y]))
                                v = i;
                } else {
                    if (rules[i].right[0].first[0].equal_symbol(non_terminal_string_map[pumping_v]))
                        if (rules[i].right[0].first[1].equal_symbol(rules[i].left[0]))
                            v = i;
                }
            } else {
                if (rules[i].right[0].first[0].equal_symbol(rules[i].left[0]))
                    if (rules[i].right[0].first[1].equal_symbol(non_terminal_string_map[pumping_y]))
                        v = i;
            }
        } else if (rules[i].right.size() >1) {
            if (rules[i].right[1].first[1].equal_symbol(non_terminal_string_map[npumping_x])) {
                if (!pumping_v.empty()) {
                    if (!pumping_y.empty()) {
                        if (rules[i].right[0].first[0].equal_symbol(non_terminal_string_map[pumping_v]))
                            if (rules[i].right[0].first[1].equal_symbol(rules[i].left[0]))
                                if (rules[i].right[0].first[2].equal_symbol(non_terminal_string_map[pumping_y]))
                                    v = i;
                    } else {
                        if (rules[i].right[0].first[0].equal_symbol(non_terminal_string_map[pumping_v]))
                            if (rules[i].right[0].first[1].equal_symbol(rules[i].left[0]))
                                v = i;
                    }
                } else {
                    if (rules[i].right[0].first[0].equal_symbol(rules[i].left[0]))
                        if (rules[i].right[0].first[1].equal_symbol(non_terminal_string_map[pumping_y]))
                            v = i;
                }
            }
        }

    }
    if (v != -1)
        return;

    if (!pumping_v.empty())
        rhs_to_look.push_back(non_terminal_string_map[pumping_v]);

    Symbol::Symbol new_nt = Symbol::Symbol("NT"+to_string(n_non_terminals), n_non_terminals, false, false);
    n_non_terminals++;
    non_terminals.push_back(new_nt);
    rhs_to_look.push_back(new_nt);



    if (!pumping_y.empty())
        rhs_to_look.push_back(non_terminal_string_map[pumping_y]);
    vector<Symbol::Symbol> lhs;
    lhs.push_back(new_nt);
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> rhs;
    rhs.push_back(make_pair(rhs_to_look, make_pair(0.0, 0.1)));

    rhs_to_look.clear();
    if (!pumping_v.empty())
        rhs_to_look.push_back(non_terminal_string_map[pumping_v]);
    if (!npumping_x.empty())
        rhs_to_look.push_back(non_terminal_string_map[npumping_x]);
    if (!pumping_y.empty())
        rhs_to_look.push_back(non_terminal_string_map[pumping_y]);
    rhs.push_back(make_pair(rhs_to_look, make_pair(0.0, 0.1)));
    rules.push_back(Rule::Rule(lhs, rhs));



    v = n_non_terminals-1;

}
std::vector<Symbol::Symbol> Grammar::Grammar::yield_string(std::vector<Symbol::Symbol> vector) {
    std::vector<Symbol::Symbol> yielded;
    std::stack<Symbol::Symbol> symbol_stack;
    for (int i = vector.size()-1; i >=0; i--)
        symbol_stack.push(vector[i]);

    while (!symbol_stack.empty()) {
        Symbol::Symbol aux = symbol_stack.top();
        symbol_stack.pop();
        if (aux.terminal)
            yielded.push_back(aux);
        else {
            for (int i = rules[aux.id].right[0].first.size()-1; i >=0; i--)
                symbol_stack.push(rules[aux.id].right[0].first[i]);
        }
    }
    return yielded;
}
void Grammar::Grammar::calculate_new_rule_from_starting_symbol(pair<vector<Symbol::Symbol>, pair<double, double>> &new_start_right, int v, std::string uvxyz, std::unordered_map<string, Symbol::Symbol> &non_terminal_string_map) {
    string u = "";
    string z = "";
    for (int i = 0; uvxyz[i] != '|'; i++)
        u.push_back(uvxyz[i]);

    for (int i = uvxyz.size()-1; uvxyz[i] != '|'; i--)
        z.push_back(uvxyz[i]);
    reverse(z.begin(),z.end());

    //std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> new_rhs;
    if (!u.empty())
        new_start_right.first.push_back(non_terminal_string_map[u]);
    new_start_right.first.push_back(rules[v].left[0]);
    if (!z.empty())
        new_start_right.first.push_back(non_terminal_string_map[z]);
    new_start_right.second.second = 0.1;
//    new_start_right.push_back(new_rhs);

}
void Grammar::Grammar::check_digram_position_integrity(unordered_map<std::string, std::tuple<int, int, int>> &digram_position) {
    for (auto m: digram_position) {
        vector <string> tokens;
        string line = m.first;
        stringstream check1(line);
        string intermediate;
        while(getline(check1, intermediate, ' '))
            tokens.push_back(intermediate);
        int flag = 1;

        int ruleOffSet = -1;
        if ( get<0>(m.second) == 0)
            ruleOffSet = -1;
        if (rules[get<0>(m.second)].right[get<2>(m.second)].first[get<1>(m.second)+ruleOffSet].name.compare(tokens[0]) != 0)
            cout << "CHECK BY M: Problem on digram: " << m.first << ". position: " << "("<<get<0>(m.second)<<","<<get<1>(m.second)<<","<<get<2>(m.second)<<")"<<endl;
        else if (rules[get<0>(m.second)].right[get<2>(m.second)].first[get<1>(m.second)+ruleOffSet+1].name.compare(tokens[1]) != 0)
            cout << "CHECK BY M: Problem on digram: " << m.first << ". position: " << "("<<get<0>(m.second)<<","<<get<1>(m.second)<<","<<get<2>(m.second)<<")"<<endl;
    }
}

void Grammar::Grammar::check_digram_position_integrity_by_rules(unordered_map<std::string, std::tuple<int, int, int>> &digram_position) {
    Rule::Rule r = rules[0];
        for(int j = 0; j < r.right.size(); j++) {
            if (r.right[j].second.first >0) {
                for (int i = 0; i < r.right[j].first.size()-1; i++) {
                    string digram =  r.right[j].first[i].name + " " +r.right[j].first[i+1].name;
                    if (digram_position.find(digram) == digram_position.end()) {
                        cout << "Digram " << digram << " not found. Should be in ("+ to_string(r.left.front().id) +","+ to_string(i+1)+","+to_string(j)+")" <<endl;
                    }

                    else if (digram.compare(rules[get<0>(digram_position[digram])].right[get<2>(digram_position[digram])].first[get<1>(digram_position[digram])-1].name + " " +
                                            rules[get<0>(digram_position[digram])].right[get<2>(digram_position[digram])].first[get<1>(digram_position[digram])].name) != 0)
                        cout << "CHECK BY R: Problem on digram: " << digram << ". position: " << "("<<get<0>(digram_position[digram])<<","<<get<1>(digram_position[digram])<<","<<get<2>(digram_position[digram])<<")"<<endl;

                }
            }

        }
}

void Grammar::Grammar::remove_unused_rules() {
    int lastI = 0;
    for (int i = 1; i < rules.size(); i++) {
        if (rules[i].right[0].second.first == 0.0) {
            lastI = i;
            int j = i+1;
            for (j = i+1; j < rules.size(); j++) {
                if (rules[j].right[0].second.first != 0.0)
                    break;
            }
            if (j == rules.size())
                break;
            rules[i].right[0].first.clear();
            rules[i].right[0].first.insert(rules[i].right[0].first.begin(), rules[j].right[0].first.begin(), rules[j].right[0].first.end());
            rules[i].right[0].second.first = rules[j].right[0].second.first;
            rules[j].right[0].second.first = 0.0;

            for (auto &r: rules) {
                for (auto &itRight: r.right) {
                    for (auto &s: itRight.first) {
                        if (s.id == j && s.terminal == false)
                            s = rules[i].left[0];
                    }
                }
            }
        }
    }
    rules.erase(rules.begin()+lastI, rules.end());
    non_terminals.erase(non_terminals.begin()+lastI, non_terminals.end());
    n_non_terminals = non_terminals.size();
}
void Grammar::Grammar::count_pumping_str_by_slice(unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, unordered_map<std::string, std::vector<std::vector<Symbol::Symbol>>> &map_pump_to_word) {

    for (int sub_amount = 0; sub_amount <word.size(); sub_amount++) {
        for (int i = 0; i <=sub_amount ; i ++) {
            vector<Symbol::Symbol> aux_word;
            aux_word.insert(aux_word.end(), word.begin()+i, word.begin()+i+word.size()-sub_amount);
            string find_pump = convert_vector_to_string(aux_word);
            Symbol::Symbol v = verify_pumping_times(Symbol::Symbol(find_pump, 0, false, false));
            int sb = 110;
            if (word.size() == 186 && sub_amount >= sb)
                cout << "Gotcha at " << sb;
            for (int sub_amount2 = i+word.size()-sub_amount; sub_amount2 <word.size(); sub_amount2++) {
                for (int j = i + word.size() - sub_amount; j <= sub_amount2; j++) {
                    vector<Symbol::Symbol> aux_word2;
                    aux_word2.insert(aux_word2.end(), word.begin() + j, word.begin() + j + word.size() - sub_amount2);
                    string find_pump2 = convert_vector_to_string(aux_word2);
                    Symbol::Symbol x = verify_pumping_times(Symbol::Symbol(find_pump2, 0, false, false));

                    if (v.name.compare("") != 0 || x.name.compare("") != 0) {
                        vector<Symbol::Symbol> u;
                        u.insert(u.end(), word.begin(), word.begin()+i);
                        string pumping = convert_vector_to_string(u);
                        u.clear();
                        u.insert(u.begin(), word.begin()+i+word.size()-sub_amount, word.begin() + j);
                        if (v.name.compare("") != 0)
                            pumping += "|" + v.name + "|" +convert_vector_to_string(u) + "|";
                        else {
                            u.insert(u.begin(), word.begin()+i, word.begin()+i+word.size()-sub_amount);
                            if (pumping.compare("") != 0)
                                pumping += " ";
                            pumping += convert_vector_to_string(u) + "|";
                        }


                        if (x.name.compare("") != 0)
                            pumping += x.name;
                        u.clear();
                        u.insert(u.begin(), word.begin() + j + word.size() - sub_amount2, word.end());
                        pumping += "|" + convert_vector_to_string(u);
                        if (v.name.compare("") == 0)
                            pumping += "||";
                        map[pumping] +=1;
                        map_pump_to_word[pumping].push_back(word);


                    }
                }
            }
        }
    }



}
std::vector<std::string> Grammar::Grammar::convert_symbol_to_vector_string(Symbol::Symbol s) {
    stringstream check1(s.name);
    string intermediate;
    vector<string> tokens;
    while(getline(check1, intermediate, ' '))
        tokens.push_back(intermediate);
    return tokens;
}
std::string Grammar::Grammar::convert_vector_to_string(std::vector<Symbol::Symbol> vs) {
    string word = "";
    if (vs.size()>0)
        word += vs[0].name;
    for (int i = 1; i <vs.size() ; i ++)
        word += " "+vs[i].name;
    return word;
}
