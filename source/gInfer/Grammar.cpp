
// Created by henrique on 02/04/19.
//NT agora é #$, NTT é #$$, p é !, pNT é !#$, #$I é #$%, I é %
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
        if (itr->freq() > 0.0) {
            std::cout << "\t";
            (*itr).print_rule();
        }
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
        non_terminals.emplace_back("#$" + std::to_string(i), i, false, false);
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

void Grammar::Grammar::train(training_method algorithm, int iterations, double alpha_alergia, double p_ratio, double time_limite) {
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
                    alergia(alpha_alergia);
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
                case pcfg_pumping_inference:
                    super_duper_pumping_inference(alpha_alergia, p_ratio, time_limite);
                    break;
                default:
                    inside_outside(iterations);
                    break;
            }
            break;
        case pcsg:
            switch(algorithm) {
                case pcsg_metropolis_hastings:
                    metropolis_hastings_pcsg(iterations, time_limite);
                    break;
                case pcsg_gibbs_sampling:
                    gibbs_sampling_pcsg(iterations);
                    break;
                default:
                    metropolis_hastings_pcsg(iterations, time_limite);
                    break;
            }
            break;
    }
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
}

void Grammar::Grammar::metropolis_hastings_pcsg(int iterations, double time_limit) {
    std::chrono::duration<double> totalTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> insideTTime = std::chrono::steady_clock::now() -  std::chrono::steady_clock::now();
    std::chrono::duration<double> insideTTimeMet = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> iterationTime = std::chrono::steady_clock::now() -  std::chrono::steady_clock::now();
    std::chrono::duration<double> newThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> updateThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> perplexityTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> remainingTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    //srand((unsigned) time(nullptr));
    size_t sentencesLength = words.size();
    auto startIt = std::chrono::steady_clock::now();
    for (size_t i = 0; i< sentencesLength; i++) {
        /*if (sentencesLength >= 10)
            if (i%(sentencesLength/10) ==0)
                std::cout << 100*(i/(1.0*sentencesLength))<< "% of trees parsed" << std::endl;*/

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
        //iterationTime += std::chrono::system_clock::now() - startIt;
    }

    int countAcceptedTress = 0;
    iterationTime = std::chrono::system_clock::now() - std::chrono::system_clock::now();
    for (int j = 0; j < iterations; j++) {
        std::mt19937 mt(rd());
        std::uniform_int_distribution<int> dist(0, sentencesLength-1);
        int p = dist(mt);
        int i = p;
        double ****insideTableKL;
        //double ****iTableN;
        actual_production.clear();
        actual_production.push_back(non_terminals[0]);
        /*if (iterations >= 5) {
            if (j%(iterations/10) == 0) {
                //auto startPerplexityTime = std::chrono::system_clock::now();
                std::pair<double,  double> pMpP = perplexity_kl(words);
                //perplexityTime += std::chrono::system_clock::now() - startPerplexityTime;
                std::cout << "   iteration: " << j << " Tree " << i <<" - Perplexity: "<< pMpP.first << " - PerplexityN: "<< pMpP.second <<" Accepted Tress: " << countAcceptedTress << " PTime: "<< perplexityTime.count() <<std::endl;
                //std::cout << "   iteration: " << j << " Tree " << i << " Accepted Tress: " << countAcceptedTress << std::endl;
            }
        }*/
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
        if (time_limit < std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startIt).count()/1000) {
            cout << "   Time Limit Reached" << endl;
            break;
        }
    }

    iterationTime = std::chrono::steady_clock::now() - startIt;
    cout << "RunTime: " << std::chrono::duration_cast<std::chrono::milliseconds>(iterationTime).count();
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
    auto newThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto updateThetaTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    auto perplexityTime = std::chrono::system_clock::now() -  std::chrono::system_clock::now();
    std::chrono::duration<double> iterationTime = std::chrono::steady_clock::now() -  std::chrono::steady_clock::now();
    auto startIt = std::chrono::steady_clock::now();
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
    iterationTime = std::chrono::steady_clock::now() - startIt;
    cout << "RunTime: " << std::chrono::duration_cast<std::chrono::milliseconds>(iterationTime).count();
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
    Symbol::Symbol nI = Symbol::Symbol("#$%", 0, false, false);
    non_terminals.push_back(nI);
    int nNT = 1;
    vector<Symbol::Symbol> vr;
    vr.push_back(nI);
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> empty_right;
    empty_right.push_back(std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>());
    empty_right[0].first.push_back(Symbol::Symbol("", -1, true, false));
    Rule::Rule r = Rule::Rule(vr, empty_right);

    rules.push_back(r);
    for (const auto& w: words) {
        string s;
        for (auto & i : w) {
            s += i.name;
            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
            bool existRule = false;
            for (itRight = rules[nI.id].right.begin(); itRight != rules[nI.id].right.end()-1; itRight++) {
                if ((*itRight).first[1].name == "#$"+s) {
                    (*itRight).second.first += 1.0;
                    existRule = true;
                    nI = (*itRight).first[1];
                    break;
                }
            }
            if (!existRule) {
                Symbol::Symbol nt = Symbol::Symbol("#$"+s, nNT, false, false);
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
                r = Rule::Rule(vr, empty_right);
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
            Symbol::Symbol nt = Symbol::Symbol("", -1, true, false);
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
    std::chrono::duration<double> iterationTime = std::chrono::steady_clock::now() -  std::chrono::steady_clock::now();
    auto startIt = std::chrono::steady_clock::now();
    non_terminals.clear();
    rules.clear();
    n_non_terminals = 0;
    gen_fpta();
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
            if (compatible_alergia(cRed, (*itBlue), alpha, rules)) {
                //cout << "merge " << cRed.name << " and " << (*itBlue).name << endl;
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

    }
    //nullify_unreacheble_rules();
    remove_unused_rules();
    for (auto & r: rules) {
        remove_unused_rule_zero_righties(r.right);
        /*for (int j = 0; j < r.right.size()-1; j ++) {
            if (r.right[j].first[0].name.empty()) {
                auto aux = r.right[r.right.size()-1];
                r.right[r.right.size()-1] = r.right[j];
                r.right[j] = aux;
            }
        }*/
    }
    generate_nt_for_t();
    normalize_probs();
    iterationTime = std::chrono::steady_clock::now() - startIt;
    cout << "RunTime: " << std::chrono::duration_cast<std::chrono::milliseconds>(iterationTime).count();
}

bool Grammar::Grammar::compatible_alergia(const Symbol::Symbol &a, const Symbol::Symbol &b, double alpha, std::vector<Rule::Rule> & vector_rules ) {
    bool correct = true;

    //nao tem os 2 ifs no alergia tradicional. coloquei para diferenciar NTs que são finais dos que não são
    /*if ((*(vector_rules[a.id].right.end() - 1)).second.first == 0 && (*(vector_rules[b.id].right.end() - 1)).second.first != 0)
        correct = false;
    if ((*(vector_rules[a.id].right.end() - 1)).second.first != 0 && (*(vector_rules[b.id].right.end() - 1)).second.first == 0)
        correct = false;*/

    if (!test_alergia((*(vector_rules[a.id].right.end() - 1)).second.first, vector_rules[a.id].freq(),
                      (*(vector_rules[b.id].right.end() - 1)).second.first, vector_rules[b.id].freq(), alpha))
        correct = false;
    for (const auto& t: terminals) {
        double dFreqa = 0.0;
        double dFreqb = 0.0;

        for (auto nt: vector_rules[a.id].right) {
            if (nt.first[0].id == t.id) {
                dFreqa =nt.second.first;
            }
        }
        for (auto nt: vector_rules[b.id].right) {
            if (nt.first[0].id == t.id) {
                dFreqb =nt.second.first;
            }
        }

        //nao tem os 2 ifs no alergia tradicional. coloquei para diferenciar NTs que possuem transição par aum terminal t dos que não possuem
        /*if (dFreqa == 0 && dFreqb != 0)
            correct = false;
        if (dFreqa != 0 && dFreqb == 0)
            correct = false;*/

        if (!test_alergia(dFreqa, vector_rules[a.id].freq(), dFreqb, vector_rules[b.id].freq(), alpha))
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
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end()-1; itRight++) {
            if ((*itRight).first[1].id == b.id) {
                stochastic_fold(a, b);
                std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> newRight;
                newRight.first.push_back((*itRight).first[0]);
                newRight.first.push_back(a);
                newRight.second.first = (*itRight).second.first;
                (*itRight).second.first = 0;
                if ((*itRule).right[(*itRule).right.size()-1].first[0].name.empty()) {
                    auto aux = (*itRule).right[(*itRule).right.size()-1];
                    (*itRule).right[(*itRule).right.size()-1] = newRight;
                    (*itRule).right.push_back(aux);
                } else
                    (*itRule).right.push_back(newRight);
                return;
            }
        }
    }
}


void Grammar::Grammar::stochastic_fold(const Symbol::Symbol& a, const Symbol::Symbol& b) {
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    for (int i = 0; i < rules[b.id].right.size()-1; i++) {


    //for (itRight = rules[b.id].right.begin(); itRight != rules[b.id].right.end()-1; itRight++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRighta;
        for (itRighta = rules[a.id].right.begin(); itRighta != rules[a.id].right.end(); itRighta++) {
            if (itRighta->first[0].id == rules[b.id].right[i].first[0].id) {
                if (itRighta->second.first >=1.0)
                    break;
            }
        }
        if (itRighta != rules[a.id].right.end()) {
            if (itRighta->second.first >= 1.0) {
                if (!itRighta->first[0].name.empty())
                    stochastic_fold(itRighta->first[1], rules[b.id].right[i].first[1]);
                (*itRighta).second.first += (rules[b.id].right[i]).second.first;
                (rules[b.id].right[i]).second.first = 0.0;
            }
        }
        else {
            if (rules[a.id].right[rules[a.id].right.size()-1].first[0].name.empty()) {
                std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> aux = rules[a.id].right[rules[a.id].right.size()-1];
                rules[a.id].right[rules[a.id].right.size()-1] = (rules[b.id].right[i]);
                rules[a.id].right.push_back(aux);
            } else
                rules[a.id].right.push_back(rules[b.id].right[i]);



            (rules[b.id].right[i]).second.first = 0.0;
        }
    //}
    }
    (*(rules[a.id].right.end()-1)).second.first += (*(rules[b.id].right.end()-1)).second.first;
    (*(rules[b.id].right.end()-1)).second.first = 0.0;
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
                if ((*itRight).first.size() > 1) {
                    notUsed.push_back((*itRight).first[1]);
                    Symbol::Symbol aux = (*itRight).first[1];
                }
                //recursive_insert_unused(notUsed, aux);

            }
            else {
                bool isUsed = false;
                for (const auto& u: used) {
                    if ((*itRight).first.size() > 1) {
                        if (itRight->first[1].id == u.id) {
                            isUsed = true;
                            break;
                        }
                    }
                }
                if (!isUsed)
                    if ((*itRight).first.size() > 1) {
                        used.push_back((*itRight).first[1]);
                    }
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
        if (freq > 0.0)
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
    Symbol::Symbol nt = Symbol::Symbol("#$%", nIds, false, false);
    non_terminals.push_back(nt);
    nt.name = "#$";
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
        Symbol::Symbol finalLHS = Symbol::Symbol("#$%", 0, false, false);
        for(unsigned long i = 0; i < w.size(); i++) {
            if (i == (unsigned long) n_non_terminals ) break;
            Symbol::Symbol nt = Symbol::Symbol("#$%", 0, false, false);
            if (i >=1)
                nt.name = "#$";
            finalLHS.name = "#$";
            for (const auto& s: nGram) {
                nt.name += s.name;
                finalLHS.name += s.name;
            }
            add_n_gram_rule_frequency(nt, w[i]);
            finalLHS.name += w[i].name;
            nGram.push_back(w[i]);
        }
        for(size_t i = n_non_terminals; i < w.size(); i++) {
            Symbol::Symbol nt = Symbol::Symbol("#$", 0, false, false);
            for (const auto& s: nGram)
                nt.name += s.name;
            add_n_gram_rule_frequency(nt, w[i]);
            finalLHS.name = nt.name;
            nGram.erase(nGram.begin());
            nGram.push_back(w[i]);
        }
        Symbol::Symbol final = Symbol::Symbol("", -1, true, false);
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
    Symbol::Symbol start = Symbol::Symbol("#$"+ to_string(n_non_terminals-1), n_non_terminals-1, false, false);
    non_terminals.push_back(start);
    vector<Symbol::Symbol> lhs;
    lhs.push_back(start);
    std::unordered_map<string, int> digram_map;
    std::unordered_map<string, tuple<int, int, int>> digram_position;
    std::unordered_map<int, int> rule_map;
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> rhss;
    rules.push_back(Rule::Rule(lhs,rhss));
    for (int i = 0; i < words.size(); i++) {
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
            if (!(dicard <1)) {
                rules[0].right[i].first.push_back(s);
                string digram = rules[0].right[i].first[rules[0].right[i].first.size() - 2].name + " " +rules[0].right[i].first[rules[0].right[i].first.size() - 1].name;

                while (verify_duplicate_digram(digram_position, digram)) {

                    if (i ==3527 || i == 3529) {
                        for (auto m: digram_position) {
                            if (get<2>(m.second) ==3527)
                                cout << "Display m: "<< m.first << ": " << get<1>(m.second) << endl;
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
    Symbol::Symbol nnt = Symbol::Symbol("#$" + to_string(n_non_terminals-1), n_non_terminals - 1, false, false);
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
    Symbol::Symbol start = Symbol::Symbol("#$"+ to_string(n_non_terminals-1), n_non_terminals-1, false, false);
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
    Symbol::Symbol nnt = Symbol::Symbol("#$" + to_string(n_non_terminals-1), n_non_terminals - 1, false, false);
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
    group_equal_rhs(rules[0].right);

    //TERMINAL RULES
    vector<Rule::Rule> rulesTerm;
    for (auto t: terminals) {
        n_non_terminals++;
        Symbol::Symbol nnt = Symbol::Symbol("#$" + to_string(n_non_terminals-1), n_non_terminals - 1, false, false);
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
                    //rulesTerm[nt_index - rules.size()].right[0].second.first += 1.0;
                    //if (itRule->left[0].id == 0)
                        rulesTerm[nt_index - rules.size()].right[0].second.first += (*itRight).second.first;
                } else {
                    vector<Symbol::Symbol> lTerminal;
                    n_non_terminals++;
                    Symbol::Symbol nnt = Symbol::Symbol("#$" + to_string(n_non_terminals-1), n_non_terminals - 1, false, false);
                    non_terminals.push_back(nnt);
                    lTerminal.push_back(nnt);
                    //rhsTerminal.push_back(make_pair(rhs_cnf, make_pair(1.0,0.1)));
                    //if (itRule->left[0].id == 0)
                        rhsTerminal.push_back(make_pair(rhs_cnf, make_pair((*itRight).second.first,0.1)));

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
                    for (auto & r: rules[iRule].right[0].first)
                        rules[r.id].right[0].second.first += (*itRight).second.first;
                    if (rules[iRule].right[0].second.first <= 0.0)
                        for (auto & r: rules[iRule].right[0].first)
                            rules[r.id].right[0].second.first -= 1;

                }
            }
        }
    }
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
    Symbol::Symbol word_0_pumped = verify_pumping_times(word[0], 1);
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
    Symbol::Symbol word_0_pumped = verify_pumping_times(word[0], 1);
    Symbol::Symbol  word_1_pumped = verify_pumping_times(word[1], 1);
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
    Symbol::Symbol word_0_pumped = verify_pumping_times(word[0], 1);
    Symbol::Symbol word_1_pumped = verify_pumping_times(word[1], 1);
    Symbol::Symbol word_2_pumped = verify_pumping_times(word[2], 1);
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
    Symbol::Symbol word_0_pumped = verify_pumping_times(word[0], 1);
    Symbol::Symbol word_1_pumped = verify_pumping_times(word[1], 1);
    Symbol::Symbol word_2_pumped = verify_pumping_times(word[2], 1);
    Symbol::Symbol word_3_pumped = verify_pumping_times(word[3], 1);
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
    Symbol::Symbol word_1_pumped = verify_pumping_times(word[1], 1);
    Symbol::Symbol word_3_pumped = verify_pumping_times(word[3], 1);
    if (word_1_pumped.name.compare("")  != 0 && word_3_pumped.name.compare("")  != 0) {
        string pumping = word[0].name +"|" +word_1_pumped.name + "|"+word[2].name +"|" + word_3_pumped.name +"|"+word[4].name;
        map[pumping] +=1;
        map_pump_to_word[pumping].push_back(original_word);
    }

}
Symbol::Symbol Grammar::Grammar::verify_pumping_times(Symbol::Symbol s, int pumping_size) {
    stringstream check1(s.name);
    string intermediate;
    vector<string> tokens;
    while(getline(check1, intermediate, ' '))
            tokens.push_back(intermediate);
    bool flag = true;
    if(tokens.size()<2)
        flag = false;
    if (tokens.size()%pumping_size == 0) {
        flag = true;
        int j = 0;
        for (j = 0; j < tokens.size()-pumping_size; j+=pumping_size) {
            for (int k = j; k <j+pumping_size; k++) {
                if (tokens[k].compare(tokens[k+pumping_size]) != 0) {
                    flag = false;
                    break;
                }
            }
            if(!flag)
                break;
        }

        if (flag == true) {
            if (j < pumping_size ) {
                s.name = "";
                return s;
            }
            s.name = "";
            for (int k = 0; k < pumping_size; k++)
                s.name += tokens[k] + " ";
            s.name = s.name.substr(0, s.name.size()-1);
            return s;
        }
    }
    else
        flag = false;
    if (flag == false)
        s.name = "";
    return s;
}

bool sortbysec(const pair<string ,int> &a,
               const pair<string,int> &b)
{
    return (a.second > b.second);
}

bool sort_p_int_double_bysec(const pair<int ,double> &a, const pair<int ,double> &b)
{
    return (a.second > b.second);
}


void Grammar::Grammar::pumping_inference(unordered_map<std::string, int> &map, unordered_map<string, vector<int>> &map_pump_to_word) {
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
        //cout << s.first << ": " << s.second << endl;
    }


    unordered_map<string, Symbol::Symbol> non_terminal_string_map;
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
        if (u.second < 1)
            break;
        vector<pair<vector<Symbol::Symbol>,pair<double, double>>>::iterator itRight;
        int v;
        /*if (countRules >= 391 && u.second == 3)
            cout << "testar Aqui. nnt = : " << n_non_terminals << endl;
        if (n_non_terminals >= 1496)
            cout << "testar Aqui. nnt = : " << n_non_terminals << endl;*/
        find_pumping_rules(v, u.first, non_terminal_string_map);
        for (itRight = rules[0].right.begin(); itRight != rules[0].right.end(); itRight++) {
            if (pumped_rule[itRight-rules[0].right.begin()] == false) {
                if (compare_pumping_use(itRight->first, (map_pump_to_word[u.first]))) {
                    countRules++;
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
                    if (itRight->second.first < 0.0)
                        cout << "negative rule itright." << endl;
                    for (auto r: itRight->first)
                        if (r.id != v) {
                            rules[r.id].right[0].second.first -= new_rhs.second.first;
                            if (rules[r.id].right[0].second.first < 0.0)
                                cout << "negative rule at the decreasing. rule: " << r.id << endl;
                        }


                    //teste de quantas regras ficam se excluirmos todas as regras por regras bombeadas.
                    //itRight->second.first = 0.0;
                    /*for (auto r: new_rhs.first) {
                        if (r.id != v)
                            rules[r.id].right[0].second.first += new_rhs.second.first;
                    }*/




                }
            }
        }
    }
    for (auto & right: new_start_right ) {
        for (auto & p: right.first) {
            if (rules[p.id].right.size() >1 ) {
                size_t i = rules[p.id].right[1].first[1].id;
                rules[i].right[0].second.first +=1.0;
            }
        }
    }
    rules[0].right.insert(rules[0].right.end(), new_start_right.begin(), new_start_right.end());
    remove_unused_rules();
    remove_unused_rule_zero_righties(rules[0].right);
    group_equal_initial_rules();
    normalize_probs();
    int count = 0;
    for (auto pr: pumped_rule) {
        if (pr == true)
            count ++;
    }
    cout << "Trues: " << count << endl;

    //para cada bombeamento com contador > que criterio (cobertura de palavras, mínimo de variáveis, %de bombeadas em relação ao tamanho do conjunto de treinamento)
        //para cada regra em NT0
            //criar nova regra do bombeamento
            //derivar palavra e comparar com o uso do bombeamento do mapa
            //se usar bombeamento
                //incrementar contador na nova regra
                //decrementar uso da regra que não será mais usada
                //marcar regra se ela usar o bombeamento.


}
bool Grammar::Grammar::compare_pumping_use(std::vector<Symbol::Symbol> vector, std::vector<int> & pumping_string) {
    std::vector<Symbol::Symbol> word = yield_string(vector);
    for (auto v: pumping_string) {
        if (equal_word(word, words[v]))
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
            pumping_v += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, next_bar-bar+1);
            bar = next_bar+1;
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
            npumping_x += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, next_bar-bar+1);
            bar = next_bar+1;
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
            pumping_y += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, next_bar-bar+1);
            bar = next_bar+1;
            next_bar = uvxyz.substr(init_pipe, next_pipe-init_pipe).find(" ", bar+1);
        }
        else {
            pumping_y += uvxyz.substr(init_pipe, next_pipe-init_pipe).substr(bar, string::npos);
            bar = next_bar;
        }
    }


    if (!pumping_v.empty())
        if (non_terminal_string_map.find(pumping_v) == non_terminal_string_map.end())
            create_new_nt_rule_substring(pumping_v, non_terminal_string_map);
    if (!npumping_x.empty())
        if (non_terminal_string_map.find(npumping_x) == non_terminal_string_map.end())
            create_new_nt_rule_substring(npumping_x, non_terminal_string_map);
    if (!pumping_y.empty())
        if (non_terminal_string_map.find(pumping_y) == non_terminal_string_map.end())
            create_new_nt_rule_substring(pumping_y, non_terminal_string_map);

    for (int i = 0; i < rules.size(); i++) {
        if (npumping_x.empty()) {
            if (!pumping_v.empty()) {
                if (!pumping_y.empty()) {
                    if (rules[i].right[0].first[0].equal_symbol(non_terminal_string_map[pumping_v]))
                        if (rules[i].right[0].first[1].equal_symbol(rules[i].left[0]))
                            if (rules[i].right[0].first[2].equal_symbol(non_terminal_string_map[pumping_y]))
                                if (rules[i].right.size() >1)
                                    if (rules[i].right[1].first.size() == 2)
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
                                    if (rules[i].right[1].first.size() == 3)
                                        v = i;
                    } else {
                        if (rules[i].right[0].first.size() == 2)
                            if (rules[i].right[0].first[0].equal_symbol(non_terminal_string_map[pumping_v]))
                                if (rules[i].right[0].first[1].equal_symbol(rules[i].left[0]))
                                    v = i;
                    }
                } else {
                    if (rules[i].right[0].first.size() == 2)
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

    Symbol::Symbol new_nt = Symbol::Symbol("#$"+to_string(n_non_terminals), n_non_terminals, false, false);
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

    if (!u.empty()) {
        if (non_terminal_string_map.find(u) == non_terminal_string_map.end())
            create_new_nt_rule_substring(u, non_terminal_string_map);
        new_start_right.first.push_back(non_terminal_string_map[u]);
        rules[non_terminal_string_map[u].id].right[0].second.first +=1.0;
    }
    new_start_right.first.push_back(rules[v].left[0]);
    if (!z.empty()) {
        if (non_terminal_string_map.find(z) == non_terminal_string_map.end())
            create_new_nt_rule_substring(z, non_terminal_string_map);
        new_start_right.first.push_back(non_terminal_string_map[z]);
        rules[non_terminal_string_map[z].id].right[0].second.first +=1.0;
    }

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
    int lastI = rules.size();
    for (int i = 1; i < rules.size(); i++) {
        if (rules[i].freq() <= 0.0) {
            lastI = i;
            int j = i+1;
            for (j = i+1; j < rules.size(); j++) {
                if (rules[j].freq() > 0.0)
                    break;
            }
            if (j == rules.size())
                break;
            rules[i].right = rules[j].right;
            rules[i].left[0].name = rules[j].left[0].name;
            for (auto &r: rules[j].right)
                r.second.first = 0.0;


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

void Grammar::Grammar::count_pumping_str_by_slice(std::unordered_map<std::string, int> &map, std::vector<Symbol::Symbol> word, unordered_map<string, vector<int>> &map_pump_to_word, int w_index, double reduction_share) {
    string aux = "";
    for (int sub_amount = 0; sub_amount <word.size(); sub_amount++) {
        for (int i = 0; i <=sub_amount ; i ++) {
            vector<Symbol::Symbol> aux_word;
            aux_word.insert(aux_word.end(), word.begin()+i, word.begin()+i+word.size()-sub_amount);
            string find_pump = convert_vector_to_string(aux_word);
            for (int pv_size = 1; pv_size < aux_word.size()/2+1; pv_size++) {
                Symbol::Symbol v = verify_pumping_times(Symbol::Symbol(find_pump, 0, false, false), pv_size);
                int v_pumping_times = aux_word.size()/pv_size;
                if (v.name.compare("") != 0 ) {
                    for (int sub_amount2 = i + word.size() - sub_amount; sub_amount2 < word.size(); sub_amount2++) {
                        for (int j = i + word.size() - sub_amount; j <= sub_amount2; j++) {

                            vector<Symbol::Symbol> u;
                            u.insert(u.end(), word.begin(), word.begin() + i);
                            aux = convert_vector_to_string(u);
                            string pumping = convert_vector_to_string(u);
                            vector<Symbol::Symbol> w;
                            //u.clear();
                            w.insert(w.begin(), word.begin() + i + word.size() - sub_amount, word.begin() + j);
                            aux = convert_vector_to_string(w);
                            pumping += "|" + v.name + "|" + convert_vector_to_string(w);
                            vector<Symbol::Symbol> aux_word2;
                            aux_word2.insert(aux_word2.end(), word.begin() + j, word.begin() + j + word.size() - sub_amount2);
                            string find_pump2 = convert_vector_to_string(aux_word2);
                            for (int px_size = 1; px_size < aux_word2.size()/2+1; px_size++) {
                                string final_pumping = pumping;
                                Symbol::Symbol x = verify_pumping_times(Symbol::Symbol(find_pump2, 0, false, false), px_size);
                                int w_empty = false;
                                if (w.empty())
                                    w_empty = true;
                                vector<Symbol::Symbol> z;
                                z.insert(z.begin(), word.begin() + j + word.size() - sub_amount2, word.end());
                                aux = convert_vector_to_string(z);
                                int x_pumping_times = 0;
                                if (x.name.compare("") != 0) {
                                    x_pumping_times = aux_word2.size() / px_size;
                                    if (x_pumping_times == v_pumping_times) {
                                        final_pumping += "|" + x.name + "|" + convert_vector_to_string(z);
                                    } else {
                                        if (!w_empty)
                                            final_pumping += " ";
                                        final_pumping += find_pump2;
                                        if (z.size() > 0) {
                                            final_pumping += " " +convert_vector_to_string(z) + "||";
                                        } else {
                                            final_pumping += "||";
                                        }

                                    }
                                } else {
                                    if (!w_empty)
                                        final_pumping += " ";
                                    final_pumping += find_pump2;
                                    if (z.size() > 0)
                                        final_pumping += " " + convert_vector_to_string(z);
                                    final_pumping+= "||";
                                }
                                int flag_used_pump = false;
                                for (auto it: map_pump_to_word[final_pumping] )
                                    if (it == w_index)
                                        flag_used_pump = true;
                                if (std::count(final_pumping.begin(), final_pumping.end(), '|') < 4)
                                    cout <<"error in final pump" << endl;
                                if (final_pumping.find("  ") != string::npos)
                                    cout << "deu merda na posicao: "<< final_pumping.find("  ") << endl;
                                if (!flag_used_pump) {
                                    //O bombeamento é contado com +1.0 ao invés de vpumping_times. Para voltar ao anterior, descomentar linha abaixo
                                    //map[final_pumping] += v_pumping_times;
                                    //if (v_pumping_times == x_pumping_times) {
                                        map[final_pumping] += 1.0;
                                        map_pump_to_word[final_pumping].push_back(w_index);
                                    //}
                                    if (((u.size()+ convert_symbol_to_vector_string(v).size() + w.size() + convert_symbol_to_vector_string(x).size() + z.size())*1.0)/word.size() < reduction_share)
                                        return;
                                }
                            }
                        }
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
void Grammar::Grammar::group_equal_initial_rules() {
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
}
void Grammar::Grammar::remove_unused_rule_zero_righties(std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> &right) {
    int lastI = right.size();
    //mudei o i inicial para 0 ao invés de 1. isso pode dar problema no pumping inference v1.
    for (int i = 0; i < right.size(); i++) {
        if (right[i].second.first <= 0.0) {
            lastI = i;
            int j = i+1;
            for (j = i+1; j < right.size(); j++) {
                if (right[j].second.first > 0.0)
                    break;
            }
            if (j == right.size())
                break;
            right[i] = right[j];
            right[j].second.first = 0.0;


        }
    }
    right.erase(right.begin()+lastI, right.end());
}
void Grammar::Grammar::create_new_nt_rule_substring(std::string sb, std::unordered_map<std::string, Symbol::Symbol> & non_terminal_string_map) {
    size_t next_bar = sb.find(" ", 0);
    size_t bar = 0;
    string aux = "";
    std::vector<Symbol::Symbol> rhs_to_look;
    if (sb.find("  ",0) != string::npos) {
        cout << "talvez deu merda na posicao " << sb.find("  ",0) << endl;
    }
    while (bar != string::npos) {
        if (next_bar != string::npos) {
            aux = sb.substr(bar, next_bar - bar);
            bar = next_bar + 1;
        } else {
            aux = sb.substr(bar, string::npos);
            bar = next_bar;
        }
        if (aux.find(" ",0) != string::npos) {
            cout << "talvez deu merda na posicao " << sb.find("  ",0) << endl;
        }
        next_bar = sb.find(" ", bar + 1);
        rhs_to_look.push_back(non_terminal_string_map[aux]);

    }
    Symbol::Symbol new_nt = Symbol::Symbol("#$"+to_string(n_non_terminals), n_non_terminals, false, false);
    n_non_terminals++;
    non_terminals.push_back(new_nt);
    vector<Symbol::Symbol> lhs;
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> rhs;
    lhs.push_back(new_nt);
    rhs.push_back(make_pair(rhs_to_look, make_pair(1, 0.1)));
    rules.push_back(Rule::Rule(lhs, rhs));
    non_terminal_string_map[sb] = new_nt;
}
/*Esse parser não funciona pq a gramática é ambígua. Como existe uma outraforma de derivar
 * uma parte bombeável, ele pode escolher o handle errado e não conseguir chegar no símbolo inicial.
 * Uma possibilidade seria derivar a partir da direita, porém regras que bombeiam duas substrings
 * podem não ser contempladas. Uma sugestão é fazer um parser top-down e tentar não deixar ele entrar em loop*/
double Grammar::Grammar::calculate_prob_word(std::vector<Symbol::Symbol> word) {
    double prob = 1.0;
    std::stack<Symbol::Symbol> symbol_stack;
    symbol_stack.push(Symbol::Symbol("ERROR", 0, false, false));
    int i = 0;
    Symbol::Symbol token = word[i];
    while (symbol_stack.top().name.compare("#$0") != 0 || i  != word.size()) {
        bool top_handle = false;
        for (int j = 1; j <= symbol_stack.size(); j++) {
            vector<Symbol::Symbol> handle;
            stack<Symbol::Symbol> aux = symbol_stack;
            for (int k = 0; k < j; k++) {
                handle.push_back(aux.top());
                aux.pop();
            }
            cout << "Handle: ";
            for (auto s: handle) {
                s.print_symbol();
                cout << " ";
            }
            cout << endl;
            Symbol::Symbol nt = is_handle(handle);

            if (!nt.name.empty()) {
                top_handle = true;
                for (int k = 0; k < handle.size(); k++)
                    symbol_stack.pop();
                symbol_stack.push(nt);
                for (auto right : rules[nt.id].right) {
                    reverse(handle.begin(), handle.end());
                    if (equal_word(handle, right.first)) {
                        prob *= right.second.first;
                        break;
                    }
                }
                break;
            }
        }
        if (!top_handle) {
            if (i < word.size()) {
                symbol_stack.push(word[i]);
                i++;
            } else {
                return 0.0;
            }
        }

    }
    return prob;
}
Symbol::Symbol Grammar::Grammar::is_handle(std::vector<Symbol::Symbol> sub_word) {
    for (auto r: rules) {
        for (auto right: r.right) {
            if (sub_word.size() == right.first.size()) {
                bool flag = true;
                for (int i = 0; i < sub_word.size(); i++) {
                    if (!sub_word[i].equal_symbol(right.first[right.first.size()-i-1])) {
                        flag = false;
                        break;
                    }
                }
                if (flag)
                    return r.left[0];

            }
        }
    }
    return Symbol::Symbol("",0,false, false);
}
double Grammar::Grammar::calculate_parse_tree_prob_top_down(std::vector<Symbol::Symbol> word) {
    double prob = 1.0;
    int right = 0;
    stack<pair<Symbol::Symbol, int>> symbol_rights_stack;
    symbol_rights_stack.push(make_pair(rules[0].left[0], 0));
    vector<Symbol::Symbol> yielded_word;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> parse_tree;
    vector<pair<int,int>> parse_tree_indexes;
    while (!symbol_rights_stack.empty() && symbol_rights_stack.top().second < rules[0].right.size()) {

        if (!symbol_rights_stack.top().first.terminal) {
            Symbol::Symbol to_yield = symbol_rights_stack.top().first;
            right = symbol_rights_stack.top().second;
            symbol_rights_stack.pop();
            for (int i = rules[to_yield.id].right[right].first.size()-1; i >=0; i--)
                symbol_rights_stack.push(make_pair(rules[to_yield.id].right[right].first[i], 0));

            parse_tree_indexes.push_back(make_pair(to_yield.id, right));

        } else {
            yielded_word.push_back(symbol_rights_stack.top().first);

            bool flag_yield_ok = true;
            for (int j = 0; j < yielded_word.size(); j++) {
                if (!yielded_word[j].equal_symbol(word[j])) {
                    flag_yield_ok = false;
                    while (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1 == parse_tree_indexes[parse_tree_indexes.size()-1].second) {
                        if (parse_tree_indexes.size() == 1 && parse_tree_indexes[parse_tree_indexes.size()-1].second == rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1)
                            break;
                        if (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first[0].terminal) {
                            if (!symbol_rights_stack.top().first.equal_symbol(yielded_word[yielded_word.size()-1]))
                                symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size()-1], 0));
                            yielded_word.pop_back();
                        }
                        for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                            symbol_rights_stack.pop();
                        symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0],parse_tree_indexes[parse_tree_indexes.size()-1].second));
                        parse_tree_indexes.pop_back();
                    }

                    if (!parse_tree_indexes.empty()) {
                        for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                            symbol_rights_stack.pop();
                        symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0], parse_tree_indexes[parse_tree_indexes.size()-1].second+1));
                        parse_tree_indexes.pop_back();
                    }
                    break;
                }
            }
            if (flag_yield_ok)
                symbol_rights_stack.pop();

        }
        if (symbol_rights_stack.empty()) {
            if (!equal_word(yielded_word, word)) {
                right = parse_tree_indexes[0].second+1;
                parse_tree_indexes.clear();
                yielded_word.clear();
                symbol_rights_stack.push(make_pair(rules[0].left[0], right));
            }
        }
    }

    if (parse_tree_indexes.empty())
        return 0.0;
    else {
        for (auto p: parse_tree_indexes)
            prob *= rules[p.first].right[p.second].second.first;

    }
    return prob;
}
void Grammar::Grammar::build_node_to_pump_map(std::unordered_map<std::string, std::vector<int>> &map_pump_to_word, std::unordered_map<int, std::vector<std::string>> & node_pump_map, int max_ntv, int max_ntw, int max_ntx, int max_ntz) {

    for (auto nt: non_terminals) {
        string u = nt.name.substr(2, string::npos);
        string aux;
        if (u.compare("%") == 0)
            u = "";
        else {
            aux =  "";
            for (int i = 0; i < u.size(); i++) {
                aux.push_back(u[i]);
                aux += " ";
            }
            u = aux.substr(0, aux.size()-1);
        }
        string z = u;
        reverse(z.begin(), z.end());
    }
}
void Grammar::Grammar::fpta_pumping_inference(unordered_map<std::string, int> &map, unordered_map<std::string, std::vector<int>> &map_pump_to_word, int max_ntv, int max_ntw, int max_ntx, int max_ntz) {
    /*
     * TODO Usar PegLib (PegParser, PegTL) no lugar do meu parser.
     *
     * TODO fazer toda busca da FPTA ser função do terminal (vai dar trabalho)
     *
     * */
    non_terminals.clear();
    rules.clear();
    n_non_terminals = 0;
    gen_fpta();
    //build_pumping_fronts();
    //int nt = 0;
   unordered_map<string, set<string>> map_nt_to_pumping_set = build_nt_pumping_sets(map, max_ntv, max_ntw, max_ntx, max_ntz);
    //std::map<pair<string,string>,bool> fpta_compatible_matrix = build_compatible_matrix();
    return;
    std::vector<Symbol::Symbol> pumped_nts;
    vector<Rule::Rule> static_rules = rules;
    vector<Symbol::Symbol> static_non_terminals = non_terminals;
    for (int i = 0; i < n_non_terminals; i++) {
        if ((i+1) % (n_non_terminals/5) == 0 )
            cout << "   Checking NT at " << (100.0*(i+1)/(n_non_terminals*1.0)) << "%" << endl;
        Symbol::Symbol nt1 = non_terminals[i];
        string u1 = nt1.name.substr(2, string::npos);
        string aux;
        if (u1.compare("%") == 0)
            u1 = "";
        else {
            aux =  "";
            for (int i = 0; i < u1.size(); i++) {
                aux.push_back(u1[i]);
                aux += " ";
            }
            u1 = aux.substr(0, aux.size()-1);
        }
        string z1 = "";
        reverse(z1.begin(), z1.end());
        //set<string> nt1_pumping_set = build_nt_pumping_set(nt1, map_pump_to_word, u1, max_ntv,  max_ntw,  max_ntx,  max_ntz);
        set<string> nt1_pumping_set = map_nt_to_pumping_set[nt1.name];
        for (int j = i+1; j < n_non_terminals; j++) {
            Symbol::Symbol nt2 = non_terminals[j];
            if (!nt1.equal_symbol(nt2)) {
                bool flag_subtree = true;
                string ntn1, ntn2;
                if (nt1.name.compare("#$%") != 0) {
                    if (nt1.name.size() < nt2.name.size()) {
                        ntn1 = nt1.name.substr(2, string::npos);
                        ntn2 = nt2.name. substr(2, string::npos);
                    }
                    else {
                        ntn1 = nt2.name. substr(2, string::npos);
                        ntn2 = nt1.name.substr(2, string::npos);
                    }
                    for (int k = 0; k < ntn1.size(); k++) {
                        if (ntn1[k] != ntn2[k]) {
                            flag_subtree = false;
                            break;
                        }

                    }
                }

                if (flag_subtree) {
                    string u2 = nt2.name.substr(2, string::npos);
                    string aux;
                    if (u2.compare("%") == 0)
                        u2 = "";
                    else {
                        aux =  "";
                        for (int i = 0; i < u2.size(); i++) {
                            aux.push_back(u2[i]);
                            aux += " ";
                        }
                        u2 = aux.substr(0, aux.size()-1);
                    }
                    //set<string> nt2_pumping_set = build_nt_pumping_set(nt2, map_pump_to_word, u2, max_ntv,  max_ntw,  max_ntx,  max_ntz);
                    set<string> nt2_pumping_set = map_nt_to_pumping_set[nt2.name];

                    //if (fpta_pumping_compatible(nt1_pumping_set, nt2_pumping_set, u1, z1, u2, 0.8)) {
                    bool nt_compatible = fpta_pumping_compatible_tree(nt1, nt2, 0.97, static_rules, static_non_terminals);
                    if (nt_compatible)
                        find_pumping_rule(nt1, nt2);
                    //bool nt_compatible = fpta_compatible_matrix[make_pair(nt1.name, nt2.name)];
                    if (nt_compatible && !nt_compatible) {
                        //cout << nt1.name<<" and " << nt2.name << " compatible: " << nt_compatible << endl;
                        Grammar g = Grammar({}, 1, {}, g.pcfg, make_pair(0, 0));
                        g.non_terminals.clear();
                        g.rules.clear();
                        g.terminals = terminals;
                        g.non_terminals.push_back(non_terminals[0]);
                        g.n_non_terminals = 1;
                        string nt1_end = "#$";
                        Symbol::Symbol nt1_s_end = non_terminals[0];
                        if (nt2.name.compare("#$%") == 0)
                            nt1_end = nt2.name;
                        else
                            nt1_end += nt2.name.substr(2,1);
                        int i= 1;
                        while (nt1_s_end.name.compare(nt1.name) != 0)  {
                            for (auto right: rules[nt1_s_end.id].right) {
                                if (right.first.size() <= 1) {
                                    if (nt1_end.compare("#$%") == 0)
                                        nt1_end = "#$";
                                    nt1_end += nt2.name.substr(2+i,1);
                                    break;
                                }
                                if (right.first[1].name.compare(nt1_end) == 0) {
                                    vector<Symbol::Symbol> left;
                                    left.push_back(nt1_s_end);
                                    nt1_s_end = right.first[1];

                                    right.first[1].id = g.n_non_terminals;
                                    g.non_terminals.push_back(right.first[1]);
                                    g.n_non_terminals++;

                                    vector<pair<vector<Symbol::Symbol>,pair<double, double>>> righties;
                                    righties.push_back(right);
                                    Rule::Rule r = Rule::Rule(left,righties);
                                    g.rules.push_back(r);
                                    if (nt1_end.compare("#$%") == 0)
                                        nt1_end = "#$";
                                    nt1_end += nt2.name.substr(2+i,1);
                                    break;
                                }

                            }
                            i++;
                        }
                        size_t index_nt1 = g.non_terminals[g.n_non_terminals-1].id;
                        pair<vector<Symbol::Symbol>,pair<double, double>> blocked_rhs;
                        for (auto rhs: rules[nt1.id].right)
                            if (rhs.first.size() > 1)
                                if (rhs.first[1].name.compare(nt2.name) == 0) {
                                    blocked_rhs = rhs; break; }

                        // REMOVIDO ABAIXO PARA FAZER O NOVO ADD SUBTREE OF NT1
                        /*if (nt2.name.substr(0, nt2.name.size()-1).size() >2) {
                            while (nt1_s_end.name.compare(nt2.name.substr(0, nt2.name.size()-1)) != 0)  {
                                for (auto right: rules[nt1_s_end.id].right) {
                                    if (right.first.size() <= 1) {
                                        if (nt1_end.compare("#$%") == 0)
                                            nt1_end = "#$";
                                        if (2+i <= nt2.name.size())
                                            nt1_end += nt2.name.substr(2+i,1);
                                        break;
                                    }
                                    if (right.first[1].name.compare(nt1_end) == 0) {
                                        vector<Symbol::Symbol> left;
                                        left.push_back(nt1_s_end);
                                        if (nt1_s_end.name.compare(nt1.name) == 0)
                                            blocked_rhs = right;
                                        nt1_s_end = right.first[1];

                                        right.first[1].id = g.n_non_terminals;
                                        g.non_terminals.push_back(right.first[1]);
                                        g.n_non_terminals++;

                                        vector<pair<vector<Symbol::Symbol>,pair<double, double>>> righties;
                                        righties.push_back(right);
                                        if (rules[nt1_s_end.id].right[rules[nt1_s_end.id].right.size()-1].first[0].name.empty())
                                            righties.push_back(rules[nt1_s_end.id].right[rules[nt1_s_end.id].right.size()-1]);
                                        Rule::Rule r = Rule::Rule(left,righties);
                                        g.rules.push_back(r);

                                        if (nt1_end.compare("#$%") == 0)
                                            nt1_end = "#$";
                                        nt1_end += nt2.name.substr(2+i,1);
                                        break;
                                    }
                                }
                                i++;
                            }
                        }
                        size_t index_nt2 = g.non_terminals[g.n_non_terminals-1].id;
                        g.build_pumping_righties(nt1_pumping_set, index_nt2, index_nt1);*/
                        size_t index_nt2 = add_subtrees_of_nt1_to_grammar(nt1, nt2, nt1_pumping_set, g, pumped_nts, index_nt1);
                        if (index_nt2 != -1 && !g.rules[index_nt2].right.empty()) {
                            //add_remaining_subtrees_to_grammar(nt1, blocked_rhs, g, pumped_nts);
                            for (auto t: g.terminals) {
                                Symbol::Symbol new_ntt = Symbol::Symbol("#$$"+t.name, g.n_non_terminals, false, false);
                                g.non_terminals.push_back(new_ntt);
                                g.n_non_terminals++;
                                for (auto & r: g.rules)
                                    for (auto & rhs: r.right) {
                                        rhs.second.first = 0.0;
                                        for (auto & s: rhs.first) {
                                            if (s.name.compare(t.name) == 0)
                                                s = new_ntt.clone();
                                            else if (s.name.compare(new_ntt.name) == 0)
                                                s = new_ntt.clone();
                                        }
                                    }
                                vector<Symbol::Symbol> left;
                                left.push_back(new_ntt);
                                pair<vector<Symbol::Symbol>,pair<double, double>> right;
                                right.first.push_back(t);
                                vector<pair<vector<Symbol::Symbol>,pair<double, double>>> righties;
                                righties.push_back(right);
                                Rule::Rule r = Rule::Rule(left,righties);
                                g.rules.push_back(r);
                            }
                            for (auto & w: words) {
                                //cout << convert_vector_to_string(w) << endl;
                                g.find_best_pumping_coverage(w);
                            }
                            for (auto & r: g.rules)
                                for (auto & rhs: r.right)
                                    for (auto & s: rhs.first)
                                        if (g.rules[s.id].right.size() == 1)
                                            if (g.rules[s.id].right[0].first.size() == 1)
                                                if (g.rules[s.id].right[0].first[0].terminal)
                                                    s = g.rules[s.id].right[0].first[0];
                            bool can_merge = false;
                            for (auto rhs: g.rules[index_nt2].right) {
                                if (rhs.first.size() > 2)
                                    if (rhs.second.first > 0.0) {
                                        can_merge = true;
                                        break;
                                    }
                            }
                            if (!check_id_nt_rule())
                                cout << "ERROR ID" << endl;
                            if (can_merge) {
                               if (pumping_merge_nt_3(nt1, nt2, g, pumped_nts))
                                   j--;
                               cout << ". Merged " << nt1.name << " and " << nt2.name << endl;
                            }
                            if (!check_id_nt_rule())
                                cout << "ERROR ID" << endl;
                            else
                                cout << endl;
                        }
                    }
                }
            }
        }
    }
    //pumping_alergia(0.01, pumped_nts);
}

std::set<std::string> Grammar::Grammar::build_nt_pumping_set(Symbol::Symbol nt, std::unordered_map<std::string, std::vector<int>> &map_pump_to_word, std::string u, int max_ntv, int max_ntw, int max_ntx, int max_ntz) {
    set<string> nt_pumping_set;
    std::vector<std::vector<Symbol::Symbol>> permutations;
    std::vector<Symbol::Symbol> word;
    vector<string> vs;
    vector<string> ws;
    vector<string> xs;
    vector<string> zs;
    vector<string> zps;
    int  n_symbols = count(u.begin(), u.end(), ' ')+1;
    if (u.empty())
        n_symbols = 0;
    for (int i = 0; i <= max_ntv || i <= max_ntw || i <= max_ntx || i <= max_ntz || i <= n_symbols; i++ )
        generate_permutation(permutations, terminals, i, word, 0);
    for (auto p: permutations) {
        string str = convert_vector_to_string(p);
        if (p.size() <= max_ntv)
            vs.push_back(str);
        if (p.size() <= max_ntw)
            ws.push_back(str);
        if (p.size() <= max_ntx)
            xs.push_back(str);
        if (p.size() <= max_ntz)
            zs.push_back(str);
        if (p.size() <= n_symbols)
            zps.push_back(str);
    }
    /*cout << nt.name << endl;*/
    for (auto & v: vs)
        for (auto & w: ws)
            for (auto & x: xs)
                for (auto & z: zs) {
                    for (auto & zp: zps) {
                        string pump = u+"|"+v+"|"+w+"|"+x+"|"+z;
                        if (!zp.empty()) {
                            if (!z.empty())
                                pump = u+"|"+v+"|"+w+"|"+x+"|"+z+ " " + zp;
                            else
                                pump = u+"|"+v+"|"+w+"|"+x+"|"+ zp;
                        }
                        /*cout << pump << endl;*/
                        if (map_pump_to_word.find(pump) != map_pump_to_word.end())
                            nt_pumping_set.insert(pump);

                    }

                }
    /*cout  << endl << endl;*/
    return nt_pumping_set;
}
bool Grammar::Grammar::fpta_pumping_compatible(std::set<std::string> snt1, std::set<std::string> snt2, std::string u1, std::string z1, std::string u2, double pump_percentage_tolerance) {

    double count_pump_compatible = 0.0;
    for (auto nt2: snt2) {

        string u1_aux = u2;
        string z2 = nt2;
        reverse(z2.begin(), z2.end());
        string pumping = z2.substr(z2.find("|"), string::npos);
        reverse(pumping.begin(), pumping.end());
        pumping = pumping.substr(pumping.find("|"), string::npos);
        z2 = z2.substr(0, z2.find("|"));

        string z1_aux = z2;
        bool flag_can_pump = true;
        if (!u1.empty())
            if (u2.find(u1) < u2.size()-2 && u2.size()-2 < u2.size())
                u1_aux = u2.substr(u2.find(u1)+2, string::npos);
            else
                flag_can_pump = false;
        if (!z1.empty())
            if (z2.find(z1) < z2.size()-2 && z2.size()-2 < z2.size())
                z1_aux = z2.substr(z2.find(z1)+2, string::npos);
            else
                flag_can_pump = false;
        size_t  n_symbols = count(u1_aux.begin(), u1_aux.end(), ' ')+1;
        int u1_aux_n_pump = 1;
        int z1_aux_n_pump = 1;
        if (z1_aux.empty() && !u1_aux.empty() || !z1_aux.empty() && u1_aux.empty())
            flag_can_pump = false;
        else if ((count(u1_aux.begin(), u1_aux.end(), ' ')+1) > (count(z1_aux.begin(), z1_aux.end(), ' ')+1)) {
            if ((count(u1_aux.begin(), u1_aux.end(), ' ')+1) % (count(z1_aux.begin(), z1_aux.end(), ' ')+1) != 0)
                flag_can_pump = false;
            else {
                u1_aux_n_pump = (count(u1_aux.begin(), u1_aux.end(), ' ') + 1) / (count(z1_aux.begin(), z1_aux.end(), ' ') + 1);
                n_symbols = count(z1_aux.begin(), z1_aux.end(), ' ')+1;
            }
        } else if ((count(u1_aux.begin(), u1_aux.end(), ' ')+1) < (count(z1_aux.begin(), z1_aux.end(), ' ')+1)) {
            if ((count(z1_aux.begin(), z1_aux.end(), ' ')+1) % (count(u1_aux.begin(), u1_aux.end(), ' ')+1) != 0)
                flag_can_pump = false;
            else {
                z1_aux_n_pump = (count(z1_aux.begin(), z1_aux.end(), ' ')+1) / (count(u1_aux.begin(), u1_aux.end(), ' ')+1);
                n_symbols = count(u1_aux.begin(), u1_aux.end(), ' ')+1;
            }
        }

        reverse(z1_aux.begin(), z1_aux.end());
        if (flag_can_pump) {
            vector<pair<string,string>> pumpeds;

            while(u1_aux.find(" ") != string::npos) {
                pumpeds.push_back(make_pair(u1_aux.substr(0, u1_aux.find(" ")), z1_aux.substr(0, z1_aux.find(" "))));
                for (int i = 0; i < u1_aux_n_pump; i++)
                    u1_aux = u1_aux.substr(u1_aux.find(" ")+1, string::npos);
                for (int i = 0; i < z1_aux_n_pump; i++)
                    z1_aux = z1_aux.substr(z1_aux.find(" ")+1, string::npos);
            }
            pumpeds.push_back(make_pair(u1_aux.substr(0, u1_aux.find(" ")), z1_aux.substr(0, z1_aux.find(" "))));
            for (int i = 1; i <= n_symbols; i++) {
                int flag_compatible = true;
                if (n_symbols%i == 0) {
                    for (int j = 0; j < n_symbols; j += i) {
                        string pump_v;
                        string pump_x;
                        for (int k = 0; k < i; k++) {
                            pump_v += pumpeds[j+k].first + " ";
                            pump_x += pumpeds[pumpeds.size()-1-j-k].second + " ";
                        }
                        pump_v = pump_v.substr(0, pump_v.size()-1);
                        pump_x = pump_x.substr(0, pump_x.size()-1);
                        string pump = "|"+pump_v+"||"+pump_x+"|";
                        if (snt1.find(u1+pump+z1) == snt1.end())
                            flag_compatible =  false;
                    }

                } else {
                    flag_compatible = false;
                }
                if (flag_compatible) {
                    if (snt1.find(u1+pumping+z1) != snt1.end()) {
                        count_pump_compatible += 1.0;
                        break;
                    }
                }
            }
        }
    }
    if (count_pump_compatible/snt1.size() > pump_percentage_tolerance) {
        cout << "YES: count: " << count_pump_compatible << ",  size: " <<snt1.size() << ", ratio: " << count_pump_compatible/snt1.size() << ".";
        return true;
    }

    cout << "NO: count: " << count_pump_compatible << ",  size: " <<snt1.size() << ", ratio: " << count_pump_compatible/snt1.size() << ".";
    return false;

}
void Grammar::Grammar::find_best_pumping_coverage(std::vector<Symbol::Symbol> word) {
    //cout << convert_vector_to_string(word) << endl;
    Grammar g  = Grammar({}, 1, {}, g.pcfg, make_pair(0, 0));
    g.rules = rules;
    g.non_terminals = non_terminals;
    g.terminals = terminals;

    if (word.empty()) {
        for (auto & rhs : rules[0].right)
            for (auto & s: rhs.first)
                if (s.name.empty()) {
                    rhs.second.first += 1.0;
                    break;
                }
        return;
    }
    int right = 0;
    stack<pair<Symbol::Symbol, int>> symbol_rights_stack;
    symbol_rights_stack.push(make_pair(rules[0].left[0], 0));
    vector<Symbol::Symbol> yielded_word;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> parse_tree;
    vector<pair<int,int>> parse_tree_indexes;
    while (!symbol_rights_stack.empty() && symbol_rights_stack.top().second < rules[symbol_rights_stack.top().first.id].right.size() || !symbol_rights_stack.empty() && symbol_rights_stack.top().second < rules[symbol_rights_stack.top().first.id].right.size() && symbol_rights_stack.top().first.terminal) {

        if (!symbol_rights_stack.top().first.terminal) {
            Symbol::Symbol to_yield = symbol_rights_stack.top().first;
            right = symbol_rights_stack.top().second;
            symbol_rights_stack.pop();
            //g.rules[to_yield.id].right[right].second.first += 1.0;
            for (int i = rules[to_yield.id].right[right].first.size()-1; i >=0; i--)
                symbol_rights_stack.push(make_pair(rules[to_yield.id].right[right].first[i], 0));

            parse_tree_indexes.push_back(make_pair(to_yield.id, right));

        } else {
            yielded_word.push_back(symbol_rights_stack.top().first);

            bool flag_yield_ok = true;
            vector<Symbol::Symbol> yielded_no_emptystr_word = remove_empty_substring(yielded_word);

            for (int j = 0; j < yielded_no_emptystr_word.size(); j++) {
                if (!yielded_no_emptystr_word[j].equal_symbol(word[j]) || yielded_no_emptystr_word.size() > word.size()) {
                    flag_yield_ok = false;
                    while (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1 == parse_tree_indexes[parse_tree_indexes.size()-1].second) {
                        if (parse_tree_indexes.size() == 1 && parse_tree_indexes[parse_tree_indexes.size()-1].second == rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1)
                            break;
                        if (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first[0].terminal) {
                            if (!symbol_rights_stack.top().first.equal_symbol(yielded_word[yielded_word.size()-1]))
                                symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size()-1], 0));
                            yielded_word.pop_back();
                        }
                        for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                            symbol_rights_stack.pop();
                        symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0],parse_tree_indexes[parse_tree_indexes.size()-1].second));
                        //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                        parse_tree_indexes.pop_back();
                    }
                    if (!yielded_word.empty()) {
                        if (yielded_word[yielded_word.size()-1].name.empty()) {
                            symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size() - 1], 0));
                            yielded_word.pop_back();
                        }
                    }
                    if (!parse_tree_indexes.empty()) {
                        for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                            symbol_rights_stack.pop();
                        symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0], parse_tree_indexes[parse_tree_indexes.size()-1].second+1));
                        //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                        parse_tree_indexes.pop_back();
                    }
                    break;
                }
            }
            if (flag_yield_ok)
                symbol_rights_stack.pop();

        }
        if (symbol_rights_stack.empty()) {
            if (!equal_word(remove_empty_substring(yielded_word), word)) {
                while (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1 == parse_tree_indexes[parse_tree_indexes.size()-1].second) {
                    if (parse_tree_indexes.size() == 1 && parse_tree_indexes[parse_tree_indexes.size()-1].second == rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1)
                        break;
                    if (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first[0].terminal) {
                        symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size()-1], 0));
                        yielded_word.pop_back();
                    }
                    for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                        symbol_rights_stack.pop();
                    symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0],parse_tree_indexes[parse_tree_indexes.size()-1].second));
                    //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                    parse_tree_indexes.pop_back();
                }
                if (!yielded_word.empty()) {
                    if (yielded_word[yielded_word.size() - 1].name.empty()) {
                        symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size() - 1], 0));
                        yielded_word.pop_back();
                    }
                }
                if (!parse_tree_indexes.empty()) {
                    for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                        symbol_rights_stack.pop();
                    symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0], parse_tree_indexes[parse_tree_indexes.size()-1].second+1));
                    //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                    parse_tree_indexes.pop_back();
                }
            }
        }
        /*if (symbol_rights_stack.top().first.name.empty()) {
            symbol_rights_stack.pop();
        }*/
    }

    if (!parse_tree_indexes.empty()) {
        for (auto p: parse_tree_indexes)
            rules[p.first].right[p.second].second.first += 1.0;

    }
    return;
}

Symbol::Symbol Grammar::Grammar::find_terminal_by_name(std::string name) {
    for (auto &t: terminals)
        if (t.name.compare(name) == 0)
            return t;
    return Symbol::Symbol("", -1, true, false);
}

Symbol::Symbol Grammar::Grammar::find_symbol_in_vector_by_name(std::string name, vector<Symbol::Symbol>  &vec) {
    for (auto &t: vec)
        if (t.name.compare(name) == 0)
            return t;
    return Symbol::Symbol("", -1, true, false);
}

std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> Grammar::Grammar::build_pumping_righties(set<string> nt1_pumping_set, size_t index_nt2, size_t index_nt1) {
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right;

    for (auto s: nt1_pumping_set) {
        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> r;
        std::vector<Symbol::Symbol> rhs_pumped;
        std::vector<Symbol::Symbol> rhs_not_pumped;
        size_t init_pipe = s.find("|", 0);
        size_t next_pipe = s.find("|", init_pipe+1);
        size_t next_bar = s.substr(init_pipe+1, next_pipe-init_pipe-1).find(" ", 0);
        string pumping_v = "";
        size_t bar = 0;

        if(next_pipe-init_pipe<=1)
            bar = string::npos;
        while (bar != string::npos) {
            if (next_bar != string::npos) {
                pumping_v = s.substr(init_pipe+1, next_pipe-init_pipe-1).substr(bar, next_bar-bar);
                rhs_pumped.push_back(find_terminal_by_name(pumping_v));
                bar = next_bar+1;
                next_bar = s.substr(init_pipe+1, next_pipe-init_pipe-1).find(" ", bar+1);
            }
            else {
                pumping_v = s.substr(init_pipe+1, next_pipe-init_pipe-1).substr(bar, string::npos);
                rhs_pumped.push_back(find_terminal_by_name(pumping_v));
                bar = next_bar;
            }
        }

        rhs_pumped.push_back(non_terminals[index_nt1]);


        init_pipe = next_pipe;
        next_pipe = s.find("|", init_pipe+1);
        string npumping_x = "";
        next_bar = s.substr(init_pipe+1, next_pipe-init_pipe-1).find(" ", 0);
        bar = 0;
        if(next_pipe-init_pipe<=1)
            bar = string::npos;
        while (bar != string::npos) {
            if (next_bar != string::npos) {
                npumping_x = s.substr(init_pipe+1, next_pipe-init_pipe-1).substr(bar, next_bar-bar);
                rhs_not_pumped.push_back(find_terminal_by_name(npumping_x));
                bar = next_bar+1;
                next_bar = s.substr(init_pipe+1, next_pipe-init_pipe-1).find(" ", bar+1);
            }
            else {
                npumping_x = s.substr(init_pipe+1, next_pipe-init_pipe-1).substr(bar, string::npos);
                bar = next_bar;
                rhs_not_pumped.push_back(find_terminal_by_name(npumping_x));
            }
        }
        if (npumping_x.empty())
            rhs_not_pumped.push_back(find_terminal_by_name(npumping_x));

        init_pipe = next_pipe;
        next_pipe = s.find("|", init_pipe+1);
        next_bar = s.substr(init_pipe+1, next_pipe-init_pipe-1).find(" ", 0);
        string pumping_y = "";
        bar = 0;
        if(next_pipe-init_pipe<=1)
            bar = string::npos;
        while (bar != string::npos) {
            if (next_bar != string::npos) {
                pumping_y = s.substr(init_pipe+1, next_pipe-init_pipe-1).substr(bar, next_bar-bar);
                rhs_pumped.push_back(find_terminal_by_name(pumping_y));
                bar = next_bar+1;
                next_bar = s.substr(init_pipe+1, next_pipe-init_pipe-1).find(" ", bar+1);
            }
            else {
                pumping_y = s.substr(init_pipe+1, next_pipe-init_pipe-1).substr(bar, string::npos);
                rhs_pumped.push_back(find_terminal_by_name(pumping_y));
                bar = next_bar;
            }
        }
        bool flag_exist_rhs = false;
        for (auto rhs: right)
            if (equal_rhs(rhs.first, rhs_pumped))
                flag_exist_rhs = true;
        if (!flag_exist_rhs) {
            r.first = rhs_pumped;
            right.push_back(r);
        }
        flag_exist_rhs = false;
        for (auto rhs: right)
            if (equal_rhs(rhs.first, rhs_not_pumped))
                flag_exist_rhs = true;
        if (!flag_exist_rhs) {
            r.first = rhs_not_pumped;
            right.push_back(r);
        }
    }
    return right;
    /*vector<Symbol::Symbol> left;

    left.push_back(non_terminals[index_nt2]);
    Rule::Rule r_pump = Rule::Rule(left, right);
    rules.push_back(r_pump);*/
}

std::vector<Symbol::Symbol> Grammar::Grammar::remove_empty_substring(std::vector<Symbol::Symbol> word) {
    for (int i = 0; i < word.size(); i ++) {
        if (word[i].name.empty()) {
            for (int j = i+1; j < word.size(); j++){
                word[j-1] = word[j];
            }
            word.pop_back();
        }
    }
    return word;
}
bool Grammar::Grammar::pumping_merge_nt(Symbol::Symbol nt1, Symbol::Symbol nt2, Grammar &g, std::vector<Symbol::Symbol> & pumped_nts) {
    string nt1_subs_str = nt1.name;
    bool flag_nt2_removed = true;
    if (nt1_subs_str.compare("#$%") == 0)
        nt1_subs_str = "#$";
    nt1_subs_str += nt2.name.substr(2,1);
    for (int i = 1; nt1.name.compare(nt2.name) != 0; nt1_subs_str += nt2.name.substr(2+i,1)) {
        for (int j = 0; j < rules[nt1.id].right.size(); j++ ) {
            if (rules[nt1.id].right[j].first.size() > 1)
                if (rules[nt1.id].right[j].first[1].name.compare(nt1_subs_str) == 0) {
                    if (!verify_nt_in_subtree_if_any_pumped(rules[nt1.id].right[j].first[1], pumped_nts))
                        rules[nt1.id].right[j].second.first = 0.0;
                    nt1 = find_symbol_in_vector_by_name(rules[nt1.id].right[j].first[1].name, non_terminals);
                    break;
                }
        }
    }
    if (rules[nt2.id].freq() > 0.0)
        flag_nt2_removed = false;
    Symbol::Symbol aux = find_symbol_in_vector_by_name(nt2.name.substr(0, nt2.name.size()-1), non_terminals);
    if (!aux.name.empty())
        pumped_nts.push_back(aux);
    queue<Symbol::Symbol> q_symbol;
    q_symbol.push(nt2);
    while (!q_symbol.empty()) {
        aux = q_symbol.front();
        if (!verify_nt_in_subtree_if_any_pumped(aux, pumped_nts)) {
            for (int i = 0; i < non_terminals.size(); i++) {
                if (non_terminals[i].name.compare(aux.name) == 0) {
                    non_terminals.erase(non_terminals.begin() + i);
                    break;
                }
            }
            q_symbol.pop();
            for (auto & rhs: rules[aux.id].right) {
                rhs.second.first = 0.0;
                if (rhs.first.size() > 1)
                    q_symbol.push(rhs.first[1]);
            }
        } else
            q_symbol.pop();
    }
    for (auto r: g.rules) {
        if (verify_nt_in_subtree_of_nt2(r.left[0], nt2) || (r.right.size() == 1 && r.right[0].first[0].terminal) ) {
            int nt_index = find_symbol_in_vector_by_name(r.left[0].name, non_terminals).id;
            if (nt_index != -1) {
                rules[nt_index].right.insert(rules[nt_index].right.end(), r.right.begin(), r.right.end());
                group_equal_rhs(rules[nt_index].right);
                remove_unused_rule_zero_righties(rules[nt_index].right);
                if (rules[nt_index].right.empty())
                    rules.erase(rules.begin() + nt_index);
            }
            else {
                if (r.freq() > 0.0) {
                    non_terminals.push_back(Symbol::Symbol(r.left[0].name, n_non_terminals, false, false));
                    n_non_terminals++;
                    rules.push_back(r);
                }
            }
        }
    }
    for (int i = 0; i < non_terminals.size(); i++)
        non_terminals[i].id = i;
    n_non_terminals = non_terminals.size();
    for (int i = 0; i < rules.size(); i++) {
        remove_unused_rule_zero_righties(rules[i].right);
        if (rules[i].right.empty()) {
            rules.erase(rules.begin()+i);
            i--;
        } else {
            bool flag_empty = false;
            for (int j =0; j < rules[i].right.size(); j++) {
                if (rules[i].right[j].first[0].name.empty())  {
                    flag_empty = true;
                    pair<vector<Symbol::Symbol>, pair<double,double>> aux = rules[i].right[rules[i].right.size()-1];
                    rules[i].right[rules[i].right.size()-1] = rules[i].right[j];
                    rules[i].right[j] = aux;
                    break;
                }
            }
            if (!flag_empty) {
                std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> empty_right;
                empty_right.first.push_back(Symbol::Symbol("", -1, true, false));
                rules[i].right.push_back(empty_right);
            }
        }
    }
    for (int i = 0; i < rules.size(); i ++) {
        for (int j = 0; j < rules[i].left.size(); j++)
            rules[i].left[j] = find_symbol_in_vector_by_name(rules[i].left[j].name, non_terminals);
        for (int j = 0; j < rules[i].right.size(); j++)
            for (int k = 0; k < rules[i].right[j].first.size(); k++)
                if (!rules[i].right[j].first[k].terminal)
                    rules[i].right[j].first[k] = find_symbol_in_vector_by_name(rules[i].right[j].first[k].name, non_terminals);
    }
    return flag_nt2_removed;
}
void Grammar::Grammar::group_equal_rhs(vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> &right) {
    for (int j = 0; j < right.size(); j++) {
        for (int i = j+1; i < right.size(); i++) {
            if (equal_rhs(right[j].first,  right[i].first)) {
                right[j].second.first += right[i].second.first;
                right.erase(right.begin()+i);
                i--;
            }
        }
    }
}

std::unordered_map<std::string, std::set<std::string>> Grammar::Grammar::build_nt_pumping_sets(unordered_map<std::string, int> & map, int max_ntv, int max_ntw, int max_ntx, int max_ntz) {
    std::unordered_map<std::string, std::set<std::string>> map_nt_to_pumping_set;
    for (auto & p: map) {
        size_t pos = p.first.find("|", 0);
        string s = "#$" + p.first.substr(0, pos);
        s.erase(remove(s.begin(), s.end(), ' '), s.end());
        if (s.size() == 2)
            s += "%";
        string aux = p.first.substr(pos+1, p.first.find("|", pos+1)- pos-1);
        aux.erase(remove(aux.begin(), aux.end(), ' '), aux.end());
        if (aux.size() <= max_ntv) {
            pos = p.first.find("|", pos+1);
            aux = p.first.substr(pos+1, p.first.find("|", pos+1)- pos-1);
            aux.erase(remove(aux.begin(), aux.end(), ' '), aux.end());
            if (aux.size() <= max_ntw) {
                pos = p.first.find("|", pos+1);
                aux = p.first.substr(pos+1, p.first.find("|", pos+1)- pos-1);
                aux.erase(remove(aux.begin(), aux.end(), ' '), aux.end());
                if (aux.size() <= max_ntx) {
                    pos = p.first.find("|", pos+1);
                    aux = p.first.substr(pos+1, p.first.size() - pos-1);
                    aux.erase(remove(aux.begin(), aux.end(), ' '), aux.end());
                    if (aux.size() <= max_ntz + p.first.substr(0, p.first.find("|", 0)).size())
                        map_nt_to_pumping_set[s].insert(p.first);
                }
            }
        }


    }
    return map_nt_to_pumping_set;
}
bool Grammar::Grammar::verify_nt_in_subtree_of_nt2(Symbol::Symbol nt1, Symbol::Symbol nt2) {
    if (nt1.name.compare("#$%") == 0)
        return true;
    if (nt1.name.size() > nt2.name.size())
        return false;
    for (int i = 0; i < nt1.name.size(); i++)
        if (nt1.name[i] != nt2.name[i])
            return false;
    return true;
}
bool Grammar::Grammar::verify_nt_in_subtree_if_any_pumped(Symbol::Symbol nt1, std::vector<Symbol::Symbol> pumped_nts) {
    for (auto s: pumped_nts)
        if (verify_nt_in_subtree_of_nt2(nt1, s))
            return  true;
    return false;

}
bool Grammar::Grammar::verify_nt_in_subtree_if_any_pumped_not_equal(Symbol::Symbol nt1, std::vector<Symbol::Symbol> pumped_nts) {
    for (auto s: pumped_nts) {
        if (!nt1.equal_symbol(s))
            if (verify_nt_in_subtree_of_nt2(nt1, s))
                return  true;
    }
    return false;

}

void Grammar::Grammar::add_remaining_subtrees_to_grammar(Symbol::Symbol nt, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> blocked_rhs, Grammar &g, std::vector<Symbol::Symbol> &pumped_nts) {
    for (auto rhs: rules[nt.id].right) {
        if (!equal_rhs(rhs.first, blocked_rhs.first)) {
            if (rhs.first.size() == 2) {
                if (rhs.first[0].terminal) {
                    g.rules[g.find_symbol_in_vector_by_name(nt.name, g.non_terminals).id].right.push_back(rhs);
                    queue<Symbol::Symbol> q_symbol;
                    q_symbol.push(rhs.first[1]);
                    while (!q_symbol.empty()) {
                        Symbol::Symbol aux = q_symbol.front();
                        q_symbol.pop();
                        g.rules.push_back(rules[aux.id]);
                        g.non_terminals.push_back(Symbol::Symbol(aux.name, g.n_non_terminals, aux.terminal, aux.context));
                        g.n_non_terminals++;
                        //if (!verify_nt_in_subtree_if_any_pumped(aux, pumped_nts)) {
                            for (auto s: rules[aux.id].right)
                                if (s.first.size() > 1) {
                                    if (verify_nt_in_subtree_of_nt2(aux, s.first[1]))
                                        q_symbol.push(s.first[1]);
                                }
                        //}
                    }
                }
            }
        }
    }
    for (int i = 0; i < g.rules.size(); i++) {
        g.rules[i].left[0] = find_symbol_in_vector_by_name(g.rules[i].left[0].name, g.non_terminals);
        for (int j = 0; j <  g.rules[i].right.size(); j++) {
            for (int k = 0; k < g.rules[i].right[j].first.size(); k++)
                if (!g.rules[i].right[j].first[k].terminal)
                    g.rules[i].right[j].first[k] = find_symbol_in_vector_by_name(g.rules[i].right[j].first[k].name, g.non_terminals);
        }
    }
}
bool Grammar::Grammar::fpta_pumping_compatible_tree(Symbol::Symbol nt1, Symbol::Symbol nt2, double tolerance, std::vector<Rule::Rule> &vector_rules, std::vector<Symbol::Symbol> &vector_symbol) {
    nt1 = find_symbol_in_vector_by_name(nt1.name, vector_symbol);
    nt2 = find_symbol_in_vector_by_name(nt2.name, vector_symbol);
    /*if (nt1.name.size() > nt2.name.size()) {
        Symbol::Symbol aux = nt1;
        nt1 = nt2;
        nt2 = aux;
    }*/
    //TODO não é simétrico. Ex: NT000 * NT0 é compatível, mas Nt0 * NT000 não é. Acredito que é pq o exist_subtree para quando não achamos  subarvore no nt2
    if ((*(vector_rules[nt1.id].right.end() - 1)).second.first == 0 && (*(vector_rules[nt2.id].right.end() - 1)).second.first != 0)
        return false;
    if ((*(vector_rules[nt1.id].right.end() - 1)).second.first != 0 && (*(vector_rules[nt2.id].right.end() - 1)).second.first == 0)
        return false;
    queue<Symbol::Symbol> q_nt1;
    queue<Symbol::Symbol> q_nt2;
    q_nt1.push(nt1);
    q_nt2.push(nt2);
    //cout << "Checking " << q_nt1.front().name <<" and " << q_nt2.front().name << endl;
    double count = 0;
    double n_nodes = 0;
    double n_sons = 0.0;
    double leaf_nodes = 0.0;
    while (!q_nt1.empty() && !q_nt2.empty()) {
        n_nodes += 1.0;
        //cout << "   " << q_nt1.front().name <<" and " << q_nt2.front().name << endl;
        if (compatible_alergia(q_nt1.front(), q_nt2.front(), tolerance, vector_rules))
            count += 1.0;
        else
            compatible_alergia(q_nt1.front(), q_nt2.front(), tolerance, vector_rules);
        if (vector_rules[q_nt1.front().id].right.empty())
            leaf_nodes += 1.0;
        for (auto rhs: vector_rules[q_nt1.front().id].right)
            if (rhs.first.size() == 2) {
                bool exist_both_subtree = false;
                for (auto rhs2: vector_rules[q_nt2.front().id].right)
                    if (rhs2.first.size() == 2)
                        if (rhs.first[0].equal_symbol(rhs2.first[0])) {
                            if (rhs.first[1].name.empty())
                                cout << "error here" << endl;
                            n_sons +=1.0;
                            q_nt1.push(rhs.first[1]);
                            q_nt2.push(rhs2.first[1]);
                            exist_both_subtree = true;
                            break;
                        }
                if (!exist_both_subtree) {
                    if (n_nodes > 4)
                        if (count/n_nodes > tolerance /*&& n_sons/n_nodes > 1.0*/) {
                            cout << nt1.name << " * " << nt2.name << " - n_nodes:  " << n_nodes << ", count: " << count << ", ratio: " <<count/n_nodes/* << ", average sons: " << n_sons/(n_nodes- leaf_nodes) */<< endl;
                            find_pumping_rule(nt1, nt2);
                            return true;
                        }
                    //cout << nt1.name << " * " << nt2.name << " - n_nodes:  " << n_nodes << ", count: " << count << ", ratio: " <<count/n_nodes << ", average sons: " << n_sons/(n_nodes - leaf_nodes) << endl;
                    return false;
                }
            }
        q_nt1.pop();
        q_nt2.pop();
    }
    //cout << nt1.name << " * " << nt2.name << " : " << count/n_nodes << endl;
    /*if (count/n_nodes >= tolerance)
        return true;*/
    return false;
}

bool Grammar::Grammar::pumping_merge_nt_2(Symbol::Symbol nt1, Symbol::Symbol nt2, Grammar &g, std::vector<Symbol::Symbol> & pumped_nts) {
    string nt1_subs_str = nt1.name;
    bool flag_nt2_removed = false;
    //if (nt1_subs_str.compare("#$%") == 0)
        nt1_subs_str = "#$";
    nt1_subs_str += nt2.name.substr(2,1);
    for (int i = 0; nt1.name.compare(nt2.name) != 0; nt1_subs_str += nt2.name.substr(2+i,1)) {
        for (int j = 0; j < rules[nt1.id].right.size(); j++ ) {
            if (rules[nt1.id].right[j].first.size() > 1)
                if (rules[nt1.id].right[j].first[1].name.compare(nt1_subs_str) == 0) {
                    if (!verify_nt_in_subtree_if_any_pumped(rules[nt1.id].right[j].first[1], pumped_nts)) {
                        if (rules[nt1.id].right[j].first[1].name.compare(nt2.name) == 0) {
                            if (g.rules[g.find_symbol_in_vector_by_name(nt1.name, g.non_terminals).id].freq() > 0.0) {
                                rules[nt1.id].right[j].second.first = 0.0;
                                flag_nt2_removed = true;
                            }
                            else
                                return false;
                        }
                    }
                    nt1 = find_symbol_in_vector_by_name(rules[nt1.id].right[j].first[1].name, non_terminals);
                    break;
                }
        }
        i++;
        if (2+i > nt2.name.size())
            break;
    }
    Symbol::Symbol aux;
    if (nt2.name.size() == 3)
        aux = find_symbol_in_vector_by_name("#$%", non_terminals);
    else
        aux = find_symbol_in_vector_by_name(nt2.name.substr(0, nt2.name.size()-1), non_terminals);
    if (!aux.name.empty())
        pumped_nts.push_back(aux);
    queue<Symbol::Symbol> q_symbol;
    q_symbol.push(nt2);
    while (!q_symbol.empty()) {
        aux = q_symbol.front();
        if (!verify_nt_in_subtree_if_any_pumped(aux, pumped_nts)) {
            for (int i = 0; i < non_terminals.size(); i++) {
                if (non_terminals[i].name.compare(aux.name) == 0) {
                    non_terminals.erase(non_terminals.begin() + i);
                    break;
                }
            }
            q_symbol.pop();
            for (auto & rhs: rules[aux.id].right) {
                rhs.second.first = 0.0;
                if (rhs.first.size() > 1)
                    q_symbol.push(rhs.first[1]);
            }
        } else
            q_symbol.pop();
    }
    for (auto r: g.rules) {
        if (verify_nt_in_subtree_of_nt2(r.left[0], nt2) /*|| (r.right.size() == 1 && r.right[0].first[0].terminal) */) {
            int nt_index = find_symbol_in_vector_by_name(r.left[0].name, non_terminals).id;
            if (nt_index != -1) {
                rules[nt_index].right.insert(rules[nt_index].right.end(), r.right.begin(), r.right.end());
                group_equal_rhs(rules[nt_index].right);
                remove_unused_rule_zero_righties(rules[nt_index].right);
                if (rules[nt_index].right.empty())
                    rules.erase(rules.begin() + nt_index);
            }
            else {
                if (r.freq() > 0.0) {
                    non_terminals.push_back(Symbol::Symbol(r.left[0].name, n_non_terminals, false, false));
                    n_non_terminals++;
                    rules.push_back(r);
                }
            }
        }
    }
    for (int i = 0; i < non_terminals.size(); i++)
        non_terminals[i].id = i;
    n_non_terminals = non_terminals.size();
    for (int i = 0; i < rules.size(); i++) {
        remove_unused_rule_zero_righties(rules[i].right);
        if (rules[i].right.empty()) {
            rules.erase(rules.begin()+i);
            i--;
        } else {
            bool flag_empty = false;
            for (int j =0; j < rules[i].right.size(); j++) {
                if (rules[i].right[j].first[0].name.empty())  {
                    flag_empty = true;
                    pair<vector<Symbol::Symbol>, pair<double,double>> aux = rules[i].right[rules[i].right.size()-1];
                    rules[i].right[rules[i].right.size()-1] = rules[i].right[j];
                    rules[i].right[j] = aux;
                    break;
                }
            }
            if (!flag_empty) {
                std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> empty_right;
                empty_right.first.push_back(Symbol::Symbol("", -1, true, false));
                rules[i].right.push_back(empty_right);
            }
        }
    }
    for (int i = 0; i < rules.size(); i ++) {
        for (int j = 0; j < rules[i].left.size(); j++)
            rules[i].left[j] = find_symbol_in_vector_by_name(rules[i].left[j].name, non_terminals);
        for (int j = 0; j < rules[i].right.size(); j++)
            for (int k = 0; k < rules[i].right[j].first.size(); k++)
                if (!rules[i].right[j].first[k].terminal)
                    rules[i].right[j].first[k] = find_symbol_in_vector_by_name(rules[i].right[j].first[k].name, non_terminals);
    }
    return flag_nt2_removed;
}
std::map<std::pair<std::string, std::string>, bool> Grammar::Grammar::build_compatible_matrix() {
    map<pair<string,string>,bool> fpta_compatible_matrix;
    for (int i = 0; i < n_non_terminals; i ++)
        for (int j = 0; j < n_non_terminals; j ++)
            fpta_compatible_matrix[make_pair(non_terminals[i].name, non_terminals[j].name)] = fpta_pumping_compatible_tree(non_terminals[i], non_terminals[j], 0.95, rules, non_terminals);
    return fpta_compatible_matrix;
}


std::map<std::vector<std::pair<int,int>>, std::map<int, double>>Grammar::Grammar::find_next_prefix_probabilities(std::vector<Symbol::Symbol> prefix) {
    map<vector<pair<int,int>>, map<int, double>> probs;
    for (auto t: terminals) {
        Symbol::Symbol new_ntt = Symbol::Symbol("#$_"+t.name, n_non_terminals, false, false);
        non_terminals.push_back(new_ntt);
        n_non_terminals++;
        for (auto & r: rules)
            for (auto & rhs: r.right) {
                for (auto & s: rhs.first) {
                    if (s.name.compare(t.name) == 0)
                        s = new_ntt.clone();
                    else if (s.name.compare(new_ntt.name) == 0)
                        s = new_ntt.clone();
                }
            }
        vector<Symbol::Symbol> left;
        left.push_back(new_ntt);
        pair<vector<Symbol::Symbol>,pair<double, double>> right;
        right.first.push_back(t);
        right.second.first = 1.0;
        vector<pair<vector<Symbol::Symbol>,pair<double, double>>> righties;
        righties.push_back(right);
        Rule::Rule r = Rule::Rule(left,righties);
        rules.push_back(r);
    }
    if (prefix.empty()) {
        for (auto & rhs : rules[0].right)
            for (auto & s: rhs.first)
                if (s.name.empty()) {
                    rhs.second.first += 1.0;
                    break;
                }
        return probs;
    }
    int right = 0;
    stack<pair<Symbol::Symbol, int>> symbol_rights_stack;
    symbol_rights_stack.push(make_pair(rules[0].left[0], 0));
    vector<Symbol::Symbol> yielded_word;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> parse_tree;
    vector<pair<int,int>> parse_tree_indexes;
    while (!symbol_rights_stack.empty() && symbol_rights_stack.top().second < rules[symbol_rights_stack.top().first.id].right.size() || !symbol_rights_stack.empty() && symbol_rights_stack.top().second < rules[symbol_rights_stack.top().first.id].right.size() && symbol_rights_stack.top().first.terminal) {

        if (!symbol_rights_stack.top().first.terminal) {
            Symbol::Symbol to_yield = symbol_rights_stack.top().first;
            right = symbol_rights_stack.top().second;
            symbol_rights_stack.pop();
            //g.rules[to_yield.id].right[right].second.first += 1.0;
            for (int i = rules[to_yield.id].right[right].first.size()-1; i >=0; i--)
                symbol_rights_stack.push(make_pair(rules[to_yield.id].right[right].first[i], 0));

            parse_tree_indexes.push_back(make_pair(to_yield.id, right));

        } else {
            yielded_word.push_back(symbol_rights_stack.top().first);

            bool flag_yield_ok = true;
            vector<Symbol::Symbol> yielded_no_emptystr_word = remove_empty_substring(yielded_word);

            for (int j = 0; j < yielded_no_emptystr_word.size(); j++) {
                if (!yielded_no_emptystr_word[j].equal_symbol(prefix[j]) || yielded_no_emptystr_word.size() > prefix.size()) {
                    flag_yield_ok = false;
                    while (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1 == parse_tree_indexes[parse_tree_indexes.size()-1].second) {
                        if (parse_tree_indexes.size() == 1 && parse_tree_indexes[parse_tree_indexes.size()-1].second == rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1)
                            break;
                        if (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first[0].terminal) {
                            if (!symbol_rights_stack.top().first.equal_symbol(yielded_word[yielded_word.size()-1]))
                                symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size()-1], 0));
                            yielded_word.pop_back();
                        }
                        for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                            symbol_rights_stack.pop();
                        symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0],parse_tree_indexes[parse_tree_indexes.size()-1].second));
                        //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                        parse_tree_indexes.pop_back();
                    }
                    if (!yielded_word.empty()) {
                        if (yielded_word[yielded_word.size()-1].name.empty()) {
                            symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size() - 1], 0));
                            yielded_word.pop_back();
                        }
                    }
                    if (!parse_tree_indexes.empty()) {
                        for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                            symbol_rights_stack.pop();
                        symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0], parse_tree_indexes[parse_tree_indexes.size()-1].second+1));
                        //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                        parse_tree_indexes.pop_back();
                    }
                    break;
                }
            }
            if (flag_yield_ok) {
                symbol_rights_stack.pop();
                if (equal_rhs(yielded_no_emptystr_word, prefix)) {
                    if (!symbol_rights_stack.empty()) {
                        for (auto rhs: rules[symbol_rights_stack.top().first.id].right) {
                            if (!rhs.first[0].name.empty())
                                if (!rhs.first[0].terminal)
                                    probs[parse_tree_indexes][rules[rhs.first[0].id].right[0].first[0].id] = rhs.second.first;
                                else
                                    probs[parse_tree_indexes][rhs.first[0].id] = rhs.second.first;
                        }
                    } else
                        probs[parse_tree_indexes][-1] = 1.0;
                }

            }

        }
        if (symbol_rights_stack.empty()) {
            if (!equal_word(remove_empty_substring(yielded_word), prefix)) {
                while (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1 == parse_tree_indexes[parse_tree_indexes.size()-1].second) {
                    if (parse_tree_indexes.size() == 1 && parse_tree_indexes[parse_tree_indexes.size()-1].second == rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1)
                        break;
                    if (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first[0].terminal) {
                        symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size()-1], 0));
                        yielded_word.pop_back();
                    }
                    for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                        symbol_rights_stack.pop();
                    symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0],parse_tree_indexes[parse_tree_indexes.size()-1].second));
                    //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                    parse_tree_indexes.pop_back();
                }
                if (!yielded_word.empty()) {
                    if (yielded_word[yielded_word.size() - 1].name.empty()) {
                        symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size() - 1], 0));
                        yielded_word.pop_back();
                    }
                }
                if (!parse_tree_indexes.empty()) {
                    for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                        symbol_rights_stack.pop();
                    symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0], parse_tree_indexes[parse_tree_indexes.size()-1].second+1));
                    //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                    parse_tree_indexes.pop_back();
                }
            }
        }
        if (symbol_rights_stack.size() == 1 && symbol_rights_stack.top().first.name.empty()) {
            probs[parse_tree_indexes][-1] = 1.0;
            break;
        }

    }
    return probs;
}
vector<pair<int, double>> Grammar::Grammar::find_prefix_ranking_probabilities(std::vector<Symbol::Symbol> prefix) {
    map<int, double> ranking;
    std::map<std::vector<std::pair<int,int>>, std::map<int, double>> next_s_by_tree = find_next_prefix_probabilities(prefix);
    double prob_all_t = 0.0;
    for (auto t: next_s_by_tree) {
        double prob_t = 1.0;
        double den_gen_s = 0.0;
        for (auto r: t.first)
            prob_t *= rules[r.first].right[r.second].second.first;
        prob_all_t += prob_t;
        for (auto p: t.second)
            den_gen_s += p.second;
        for (auto p: t.second)
            ranking[p.first] += (p.second/den_gen_s) * prob_t;
    }
    vector<pair<int, double>> ordered_ranking;
    for (auto w : ranking) {
        w.second /= prob_all_t;
        ordered_ranking.push_back(w);
    }
    sort(ordered_ranking.begin(), ordered_ranking.end(), sort_p_int_double_bysec);
    return ordered_ranking;
}
size_t Grammar::Grammar::add_subtrees_of_nt1_to_grammar(Symbol::Symbol nt1, Symbol::Symbol nt2, set<string> &nt1_pumping_set, Grammar &g, vector<Symbol::Symbol> &pumped_nts, size_t index_nt1) {
    size_t index_nt2 = -1;
    queue<Symbol::Symbol> q_symbol;
    q_symbol.push(nt1);
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> rhss;
    while (!q_symbol.empty()) {
        Symbol::Symbol aux = q_symbol.front();
        q_symbol.pop();
        g.rules.push_back(rules[aux.id]);
        auto rules_it = rules[aux.id].right.end();
        int i = 0;
        for (auto & s: rules[aux.id].right)
            if (s.first.size() > 1) {
                if (s.first[1].equal_symbol(nt2)) {
                    index_nt2 = find_symbol_in_vector_by_name(aux.name, g.non_terminals).id;
                    rhss = g.build_pumping_righties(nt1_pumping_set, index_nt2, index_nt1);
                    group_equal_rhs(rhss);
                    rules_it = g.rules[index_nt2].right.begin()+i;
                }
                else if (verify_nt_in_subtree_of_nt2(aux, s.first[1]) && s.first[1].name.compare("#$%") != 0) {
                    g.non_terminals.push_back(Symbol::Symbol(s.first[1].name, g.n_non_terminals, s.first[1].terminal, s.first[1].context));
                    g.n_non_terminals++;
                    q_symbol.push(s.first[1]);
                }
                i++;
            }
        if (find_symbol_in_vector_by_name(aux.name, g.non_terminals).id == index_nt2) {
            g.rules[index_nt2].right.erase(rules_it);
            rhss.insert(rhss.end(), g.rules[index_nt2].right.begin(), g.rules[index_nt2].right.end());
            group_equal_rhs(rhss);
            g.rules[index_nt2].right = rhss;
        }
    }
    for (int i = 0; i < g.rules.size(); i++) {
        g.rules[i].left[0] = find_symbol_in_vector_by_name(g.rules[i].left[0].name, g.non_terminals);
        for (int j = 0; j <  g.rules[i].right.size(); j++) {
            for (int k = 0; k < g.rules[i].right[j].first.size(); k++)
                if (!g.rules[i].right[j].first[k].terminal)
                    g.rules[i].right[j].first[k] = find_symbol_in_vector_by_name(g.rules[i].right[j].first[k].name, g.non_terminals);
        }
    }
    return index_nt2;
}

bool Grammar::Grammar::pumping_merge_nt_3(Symbol::Symbol nt1, Symbol::Symbol nt2, Grammar &g, std::vector<Symbol::Symbol> &pumped_nts) {

    Symbol::Symbol fixed_nt1 = nt1;
    string nt1_subs_str = nt1.name;
    bool flag_nt2_removed = false;
    //if (nt1_subs_str.compare("#$%") == 0)
    nt1_subs_str = "#$";
    nt1_subs_str += nt2.name.substr(2,1);
    for (int i = 0; nt1.name.compare(nt2.name) != 0; nt1_subs_str += nt2.name.substr(2+i,1)) {
        for (int j = 0; j < rules[nt1.id].right.size(); j++ ) {
            if (rules[nt1.id].right[j].first.size() == 2 )
                if (rules[nt1.id].right[j].first[1].name.compare(nt1_subs_str) == 0) {
                    if (!verify_nt_in_subtree_if_any_pumped(rules[nt1.id].right[j].first[1], pumped_nts)) {
                        if (rules[nt1.id].right[j].first[1].name.compare(nt2.name) == 0) {
                            flag_nt2_removed = true;
                        }
                    }
                    fold_fpta_subtree(fixed_nt1, nt2, pumped_nts, rules, non_terminals);
                    if (flag_nt2_removed)
                        rules[nt1.id].right[j].second.first = 0.0;
                    nt1 = find_symbol_in_vector_by_name(rules[nt1.id].right[j].first[1].name, non_terminals);
                    break;
                }
        }
        i++;
        if (2+i > nt2.name.size())
            break;
    }
    Symbol::Symbol aux;
    if (nt2.name.size() == 3)
        aux = find_symbol_in_vector_by_name("#$%", non_terminals);
    else
        aux = find_symbol_in_vector_by_name(nt2.name.substr(0, nt2.name.size()-1), non_terminals);
    if (!aux.name.empty())
        pumped_nts.push_back(aux);
    queue<Symbol::Symbol> q_symbol;
    q_symbol.push(nt2);
    while (!q_symbol.empty()) {
        aux = q_symbol.front();
        if (!verify_nt_in_subtree_if_any_pumped(aux, pumped_nts)) {
            for (int i = 0; i < non_terminals.size(); i++) {
                if (non_terminals[i].name.compare(aux.name) == 0) {
                    non_terminals.erase(non_terminals.begin() + i);
                    break;
                }
            }
        }
            q_symbol.pop();
            for (auto & rhs: rules[aux.id].right) {
                if (rhs.first.size() == 1)
                    rhs.second.first = 0.0;
                else {
                    if (rhs.first.size() > 1) {
                        if (verify_nt_in_subtree_of_nt2(aux, rhs.first[1]) && rhs.first[1].name.compare("#$%") != 0) {
                            rhs.second.first = 0.0;
                            q_symbol.push(rhs.first[1]);
                        }
                    }
                }
            }
    }
    for (auto r: g.rules) {
        if (verify_nt_in_subtree_if_any_pumped(r.left[0], pumped_nts) /*|| (r.right.size() == 1 && r.right[0].first[0].terminal) */) {
            int nt_index = find_symbol_in_vector_by_name(r.left[0].name, non_terminals).id;
            if (nt_index != -1) {
                for (auto & rhs: r.right) {
                    bool exist_both_subtree = false;
                    for (auto & rhs2: rules[nt_index].right) {
                        if (equal_rhs(rhs.first, rhs2.first)) {
                            if (rhs.first.size() == 1 || verify_nt_in_subtree_if_any_pumped(rhs.first[1], pumped_nts)) {
                                if (rhs.first.size() != 2)
                                    rhs2.second.first = rhs.second.first;
                            }
                            exist_both_subtree = true;
                            break;
                        }
                    }
                    if (!exist_both_subtree) {
                        if (rhs.first.size() == 1 || verify_nt_in_subtree_if_any_pumped(rhs.first[1], pumped_nts)) {
                            rules[nt_index].right.push_back(rhs);
                        }
                    }
                }
            } else {
                if (r.freq() > 0.0) {
                    non_terminals.push_back(Symbol::Symbol(r.left[0].name, n_non_terminals, false, false));
                    n_non_terminals++;
                    rules.push_back(r);
                }
            }
        }
    }
    /*aux = find_symbol_in_vector_by_name(nt2.name.substr(0, nt2.name.size()-1), non_terminals);
    rules[aux.id].right.insert(rules[aux.id].right.end(),
                               g.rules[find_symbol_in_vector_by_name(nt2.name.substr(0, nt2.name.size()-1), non_terminals).id].right.begin(),
                               g.rules[find_symbol_in_vector_by_name(nt2.name.substr(0, nt2.name.size()-1), non_terminals).id].right.end());
    group_equal_rhs(rules[aux.id].right);
    remove_unused_rule_zero_righties(rules[aux.id].right);*/
    for (int i = 0; i < non_terminals.size(); i++)
        non_terminals[i].id = i;
    n_non_terminals = non_terminals.size();
    for (int i = 0; i < rules.size(); i++) {

        remove_unused_rule_zero_righties(rules[i].right);
        if (rules[i].right.empty()) {
            if (!verify_nt_in_subtree_if_any_pumped(rules[i].left[0], pumped_nts)) {
            rules.erase(rules.begin()+i);
            i--;
            }
        } else {
            bool flag_empty = false;
            for (int j =0; j < rules[i].right.size(); j++) {
                if (rules[i].right[j].first[0].name.empty())  {
                    flag_empty = true;
                    pair<vector<Symbol::Symbol>, pair<double,double>> aux = rules[i].right[rules[i].right.size()-1];
                    rules[i].right[rules[i].right.size()-1] = rules[i].right[j];
                    rules[i].right[j] = aux;
                    break;
                }
            }
            if (!flag_empty) {
                std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> empty_right;
                empty_right.first.push_back(Symbol::Symbol("", -1, true, false));
                rules[i].right.push_back(empty_right);
            }
        }
    }
    for (int i = 0; i < rules.size(); i ++) {
        for (int j = 0; j < rules[i].left.size(); j++)
            rules[i].left[j] = find_symbol_in_vector_by_name(rules[i].left[j].name, non_terminals);
        for (int j = 0; j < rules[i].right.size(); j++)
            for (int k = 0; k < rules[i].right[j].first.size(); k++)
                if (!rules[i].right[j].first[k].terminal)
                    rules[i].right[j].first[k] = find_symbol_in_vector_by_name(rules[i].right[j].first[k].name, non_terminals);
    }
    return flag_nt2_removed;
}
void Grammar::Grammar::fold_fpta_subtree(Symbol::Symbol nt1, Symbol::Symbol nt2, vector<Symbol::Symbol> &pumped_nts, vector<Rule::Rule> &vector_rules, vector<Symbol::Symbol> &vector_symbol) {
    nt1 = find_symbol_in_vector_by_name(nt1.name, vector_symbol);
    nt2 = find_symbol_in_vector_by_name(nt2.name, vector_symbol);
    queue<Symbol::Symbol> q_nt1;
    queue<Symbol::Symbol> q_nt2;
    queue<Symbol::Symbol> q_nt2_remaining;
    q_nt1.push(nt1);
    q_nt2.push(nt2);
    double count_pump = 0.0;
    while (!q_nt1.empty() && !q_nt2.empty()) {
        for (auto & rhs: vector_rules[q_nt2.front().id].right) {
            bool exist_both_subtree = false;
            for (auto & rhs2: vector_rules[q_nt1.front().id].right) {
                if (rhs.first[0].id == rhs2.first[0].id) {
                    if (rhs.first.size() == 1 || !verify_nt_in_subtree_if_any_pumped(rhs.first[1], pumped_nts)) {
                        if (rhs2.first[1].equal_symbol(nt2))
                            count_pump += rhs.second.first;
                        rhs2.second.first += rhs.second.first;
                        rhs.second.first = 0.0;
                    }
                    if (rhs.first.size() > 1) {
                        if (verify_nt_in_subtree_of_nt2(q_nt2.front(), rhs.first[1]) && rhs.first[1].name.compare("#$%") != 0) {
                            q_nt2.push(rhs.first[1]);
                            if (rhs2.first[1].equal_symbol(nt2)) {
                                q_nt1.push(nt1);
                                //rhs2.second.first = 0.0;
                            }
                            else
                                q_nt1.push(rhs2.first[1]);
                        }
                    }
                    exist_both_subtree = true;
                    break;
                }
            }
            if (!exist_both_subtree) {
                q_nt2_remaining.push(rhs.first[1]);
                //COMENTEI TUDO PQ ACREDITO QUE CASO A SUBTREE NÃO EXISTA, OS SIMBÓLIS DELA SÃO GERADOS PELO BOMBEAMENTO
                /*auto right = rhs;
                Symbol::Symbol new_nt = Symbol::Symbol(q_nt1.front().name + right.first[0].name, n_non_terminals, false, false);
                n_non_terminals++;
                right.first[1] = new_nt;
                vector_rules[q_nt1.front().id].right.push_back(right);

                vector<Symbol::Symbol> lhs;
                lhs.push_back(new_nt);
                std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> new_rhs;
                Rule::Rule r = Rule::Rule(lhs, new_rhs);*/

                if (rhs.first.size() == 1 || !verify_nt_in_subtree_if_any_pumped(rhs.first[1], pumped_nts)) {
                    rhs.second.first = 0.0;
                }
                /*if (rhs.first.size() > 1) {
                    if (verify_nt_in_subtree_of_nt2(q_nt2.front(), rhs.first[1]) && rhs.first[1].name.compare("#$%") != 0) {
                        q_nt1.push(new_nt);
                        q_nt2.push(rhs.first[1]);
                    }
                }*/
            }
        }
        q_nt1.pop();
        q_nt2.pop();
    }
    while (!q_nt2_remaining.empty()) {
        for (auto & rhs: vector_rules[q_nt2_remaining.front().id].right) {
            if (rhs.first.size() == 1 || !verify_nt_in_subtree_if_any_pumped(rhs.first[1], pumped_nts)) {
                rhs.second.first = 0.0;
            }
            if (rhs.first.size() > 1) {
                if (verify_nt_in_subtree_of_nt2(q_nt2_remaining.front(), rhs.first[1]) && rhs.first[1].name.compare("#$%") != 0) {
                    q_nt2_remaining.push(rhs.first[1]);
                }
            }
        }
        q_nt2_remaining.pop();
    }
    //cout << "coutn pump: " << count_pump << endl;
}
void Grammar::Grammar::pumping_alergia(double alpha, vector<Symbol::Symbol> & pumped_nts) {
    vector<Symbol::Symbol> red;
    red.push_back(non_terminals[0]);
    vector<Symbol::Symbol> blue;
    for (auto s: rules[0].right) {
        if (s.first[0].name.empty())
            break;
        blue.push_back(s.first[1]);
    }
    int t0 = 10;


    auto itBlue = blue.begin();
    double freqCBlue = rules[(*itBlue).id].freq();
    while (!blue.empty()) {
        if (freqCBlue >= t0) {
            bool existRed = false;
            for (const auto& cRed: red) {
                if (fpta_single_pumping_compatible_tree(cRed, (*itBlue), 0.97, rules, non_terminals)) {
                    cout <<  ". merge " << cRed.name << " and " << (*itBlue).name << endl;
                    stochastic_pumping_merge(cRed, (*itBlue), pumped_nts);
                    existRed = true;
                    break;
                }
            }
            if (!existRed) {
                red.push_back((*itBlue));
                for (auto right: rules[(*itBlue).id].right) {
                    if (!(right.first[0].name.empty()) && right.first.size() == 2)
                        blue.push_back(right.first[1]);
                }
                itBlue = blue.begin();
            }
            blue.erase(itBlue);
            itBlue = blue.begin();
            if (itBlue == blue.end())
                break;
            freqCBlue = rules[(*itBlue).id].freq();
        } else {
            itBlue = blue.begin();
            blue.erase(itBlue);
            itBlue = blue.begin();
            if (itBlue == blue.end())
                break;
            freqCBlue = rules[(*itBlue).id].freq();
        }

    }
    remove_unused_rules();
    normalize_probs();
}
void Grammar::Grammar::stochastic_pumping_merge(const Symbol::Symbol &a, const Symbol::Symbol &b, vector<Symbol::Symbol> &pumped_nts) {
    vector<Rule::Rule>::iterator itRule;
    for (itRule = rules.begin(); itRule < rules.end(); itRule++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
        for (itRight = (*itRule).right.begin(); itRight != (*itRule).right.end(); itRight++) {
            if (itRight->first.size() == 2) {
                if ((*itRight).first[1].id == b.id) {
                    double count_pumping_rule = stochastic_pumping_fold(a, b, pumped_nts);
                    std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> newRight;
                    //int n = (*itRight).second.first;
                    newRight.first.push_back((*itRight).first[0]);
                    newRight.first.push_back(a);
                    newRight.second.first = (*itRight).second.first - count_pumping_rule;
                    (*itRight).second.first = count_pumping_rule;
                    (*itRule).right.insert((*itRule).right.end()-1, newRight);
                    return;
                }
            }
        }
    }
}
double Grammar::Grammar::stochastic_pumping_fold(const Symbol::Symbol &a, const Symbol::Symbol &b, vector<Symbol::Symbol> &pumped_nts) {
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    double count_pumping_rules = 0.0;
    for (itRight = rules[b.id].right.begin(); itRight != rules[b.id].right.end()-1; itRight++) {
        std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRighta;
        for (itRighta = rules[a.id].right.begin(); itRighta != rules[a.id].right.end(); itRighta++) {
            if (itRighta->first.size() <= 2 && itRight->first.size() <= 2)
                if (itRighta->first[0].id == itRight->first[0].id) {
                    if (itRighta->second.first >=1.0) {
                        if (!verify_nt_in_subtree_if_any_pumped(itRighta->first[1], pumped_nts))
                            break;
                        else {
                            count_pumping_rules += itRighta->second.first;
                        }

                    }
                }
        }
        if (itRighta != rules[a.id].right.end()) {
            if (itRighta->second.first >= 1.0) {
                if (!itRighta->first[0].name.empty())
                    stochastic_pumping_fold(itRighta->first[1], itRight->first[1], pumped_nts);
                (*itRighta).second.first += (*itRight).second.first;
                (*itRight).second.first = 0.0;
            }
        }
        else {
            rules[a.id].right.push_back(*itRight);
            (*itRight).second.first = 0.0;
        }
    }
    (*(rules[a.id].right.end()-1)).second.first += (*(rules[b.id].right.end()-1)).second.first;
    (*(rules[b.id].right.end()-1)).second.first = 0.0;
    return count_pumping_rules;
}
bool Grammar::Grammar::check_id_nt_rule() {
    for (int i = 0; i <  rules.size();i ++)
        if (rules[i].left[0].id != i)
            return false;
    return true;
}


bool Grammar::Grammar::fpta_single_pumping_compatible_tree(Symbol::Symbol nt1, Symbol::Symbol nt2, double tolerance, std::vector<Rule::Rule> &vector_rules, std::vector<Symbol::Symbol> &vector_symbol) {
    nt1 = find_symbol_in_vector_by_name(nt1.name, vector_symbol);
    nt2 = find_symbol_in_vector_by_name(nt2.name, vector_symbol);
    queue<Symbol::Symbol> q_nt1;
    queue<Symbol::Symbol> q_nt2;
    q_nt1.push(nt1);
    q_nt2.push(nt2);
    double count = 0;
    double n_nodes = 0;
    double n_sons = 0.0;
    double leaf_nodes = 0.0;
    while (!q_nt1.empty() && !q_nt2.empty()) {
        n_nodes += 1.0;
        if (compatible_alergia(q_nt1.front(), q_nt2.front(), tolerance, vector_rules))
            count += 1.0;
        else
            compatible_alergia(q_nt1.front(), q_nt2.front(), tolerance, vector_rules);
        if (vector_rules[q_nt1.front().id].right.empty())
            leaf_nodes += 1.0;
        for (auto rhs: vector_rules[q_nt1.front().id].right)
            if (rhs.first.size() == 2) {
                bool exist_both_subtree = false;
                for (auto rhs2: vector_rules[q_nt2.front().id].right)
                    if (rhs2.first.size() == 2)
                        if (rhs.first[0].equal_symbol(rhs2.first[0])) {
                            if (rhs.first[1].name.empty())
                                cout << "error here" << endl;
                            n_sons +=1.0;
                            q_nt1.push(rhs.first[1]);
                            q_nt2.push(rhs2.first[1]);
                            exist_both_subtree = true;
                            break;
                        }
                if (!exist_both_subtree) {
                    if (n_nodes > 4)
                        if (count/n_nodes > tolerance && n_sons/n_nodes <= 1.0) {
                            cout << nt1.name << " * " << nt2.name << " - n_nodes:  " << n_nodes << ", count: " << count << ", ratio: " <<count/n_nodes << ", average sons: " << n_sons/(n_nodes- leaf_nodes);
                            return true;
                        }
                    //cout << nt1.name << " * " << nt2.name << " - n_nodes:  " << n_nodes << ", count: " << count << ", ratio: " <<count/n_nodes << ", average sons: " << n_sons/(n_nodes - leaf_nodes) << endl;
                    return false;
                }
            }
        q_nt1.pop();
        q_nt2.pop();
    }
    cout << nt1.name << " * " << nt2.name << " : " << count/n_nodes << endl;
    if (count/n_nodes >= tolerance)
        return true;
    return false;
}
void Grammar::Grammar::eliminate_covered_pumpings(unordered_map<std::string, std::vector<int>> &map_pump_to_word, unordered_map<std::string, int> &map) {
    for (auto & p1: map_pump_to_word) {
        for (auto & p2: map_pump_to_word) {
            if (p1.first.compare(p2.first) != 0) {
                set<int> p1_set (p1.second.begin(), p1.second.end());
                set<int> p2_set (p2.second.begin(), p2.second.end());
                bool contain = true;
                for (auto w: p2_set) {
                    if (p1_set.find(w) == p1_set.end()) {
                        contain = false;
                        break;
                    }
                }
                if (contain) {
                    map_pump_to_word.erase(p2.first);
                    map.erase(p2.first);
                    break;
                }
                for (auto w: p1_set) {
                    if (p2_set.find(w) == p2_set.end()) {
                        contain = false;
                        break;
                    }
                }
                if (contain) {
                    map_pump_to_word.erase(p1.first);
                    map.erase(p1.first);
                    break;
                }
            }
        }
    }
}

void Grammar::Grammar::find_pumping_rule(Symbol::Symbol nt1, Symbol::Symbol nt2) {
    for (auto r: rules) {
        for (auto rhs: r.right) {
            std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right;
            std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rhs_pumping;
            if (rhs.first.size() > 1) {
                if (rhs.first[1].id == nt2.id) {
                    rhs_pumping.first.push_back(rhs.first[0]);
                    rhs_pumping.first.push_back(nt1);
                    rhs_pumping.second.first = 1.0;
                    build_rhs_until_empty_rule(nt2, right, rhs_pumping);
                    build_rhs_until_empty_rule(nt1, right, rhs_pumping);
                    vector<Symbol::Symbol> lhs = r.left;
                    Rule::Rule new_rule = Rule::Rule(lhs, right);
                    new_rule.print_rule();

                }
            }
        }
    }
}
bool Grammar::Grammar::exist_empty_rule(std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right) {
    for (auto rhs: right)
        if (rhs.first[0].name.empty())
            if (rhs.second.first > 0.0)
                return true;
    return false;
}
void Grammar::Grammar::build_rhs_until_empty_rule(Symbol::Symbol nt1, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> &right, std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs_pumping) {
    if (!exist_empty_rule(rules[nt1.id].right)) {
        for (auto rhs: rules[nt1.id].right) {
            if (!rhs.first[0].name.empty()) {
                std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs_pumping_aux = rhs_pumping;
                rhs_pumping_aux.first.push_back(rhs.first[0]);
                build_rhs_until_empty_rule(rhs.first[1], right, rhs_pumping_aux);
            }
        }
    } else {
        right.push_back(rhs_pumping);
    }
}
void Grammar::Grammar::build_pumping_fronts() {
    vector<set<int>> pumping_fronts;

    pair<int, int> first_i_j = find_first_compatible_pump();

    while (first_i_j.first != -1) {
        set<int> pumping_front;
        pumping_front.insert(first_i_j.first);
        if (first_i_j.second != -1) {
            pumping_front.insert(first_i_j.second);
            set<int> front_aux = pumping_front;
            while (!front_aux.empty()) {
                set<int> next_nts_in_front;
                for (auto nt : front_aux) {
                    for (int k = non_terminals[nt].id+1; k < n_non_terminals; k++) {
                        if (fpta_pumping_compatible_tree(non_terminals[nt], non_terminals[k], 0.95, rules, non_terminals)) {
                            if (pumping_front.find(non_terminals[k].id) == pumping_front.end()) {
                                next_nts_in_front.insert(k);
                            }
                        }
                    }
                }
                pumping_front.insert(next_nts_in_front.begin(), next_nts_in_front.end());
                front_aux = next_nts_in_front;
            }

        }
        pumping_fronts.push_back(pumping_front);

        first_i_j.first = first_i_j.second = -1;
        for (int k = 0; k < n_non_terminals; k++) {
            bool found_k = true;
            for (auto pf : pumping_fronts) {
                if (pf.find(k) != pf.end()) {
                    found_k = false;
                    break;
                }
            }
            if (found_k) {
                first_i_j.first = k;
                break;
            }
        }
        if (first_i_j.first != -1)
            for (int i = first_i_j.first + 1; i < n_non_terminals; i++) {
                if (fpta_pumping_compatible_tree(non_terminals[first_i_j.first], non_terminals[i], 0.95, rules, non_terminals)) {
                    first_i_j.second = i;
                    break;
                }

            }
    }

    for (int i=0; i < pumping_fronts.size(); i++) {
        cout << "front " << i;
        for (auto nt: pumping_fronts[i])
            cout << " " << non_terminals[nt].name << " ";
        cout << endl;
    }
    for (auto s: pumping_fronts) {
        for (auto f:  s) {
            if (*s.begin() != f) {
                cout << non_terminals[*s.begin()].name << " " << non_terminals[f].name << endl;
                find_pumping_rule(non_terminals[*s.begin()], non_terminals[f]);
            }
        }
    }


}
std::pair<int, int> Grammar::Grammar::find_first_compatible_pump() {
    for (int i = 0; i < n_non_terminals; i++)
        for (int j = i+1; j < n_non_terminals; j++)
            if (fpta_pumping_compatible_tree(non_terminals[i], non_terminals[j], 0.95, rules, non_terminals))
                return make_pair(i,j);
    return make_pair(-1,-1);
}

bool Grammar::Grammar::fpta_pumping_compatible_tree2(Symbol::Symbol nt1, Symbol::Symbol nt2, double tolerance, std::vector<Rule::Rule> &vector_rules, std::vector<Symbol::Symbol> &vector_symbol) {
    nt1 = find_symbol_in_vector_by_name(nt1.name, vector_symbol);
    nt2 = find_symbol_in_vector_by_name(nt2.name, vector_symbol);
    /*if (nt1.name.size() > nt2.name.size()) {
        Symbol::Symbol aux = nt1;
        nt1 = nt2;
        nt2 = aux;
    }*/
    if ((*(vector_rules[nt1.id].right.end() - 1)).second.first == 0 && (*(vector_rules[nt2.id].right.end() - 1)).second.first != 0)
        return false;
    if ((*(vector_rules[nt1.id].right.end() - 1)).second.first != 0 && (*(vector_rules[nt2.id].right.end() - 1)).second.first == 0)
        return false;
    queue<Symbol::Symbol> q_nt1;
    queue<Symbol::Symbol> q_nt2;
    q_nt1.push(nt1);
    q_nt2.push(nt2);
    //cout << "Checking " << q_nt1.front().name <<" and " << q_nt2.front().name << endl;
    double count = 0;
    double n_nodes = 0;
    double n_sons = 0.0;
    double leaf_nodes = 0.0;
    while (!q_nt1.empty() && !q_nt2.empty()) {
        n_nodes += 1.0;
        //cout << "   " << q_nt1.front().name <<" and " << q_nt2.front().name << endl;
        if (compatible_alergia(q_nt1.front(), q_nt2.front(), tolerance, vector_rules))
        //if (compatible_alergia(q_nt1.front(), q_nt2.front(), tolerance, vector_rules))
            count += 1.0;
        else
            return false;
        if (vector_rules[q_nt1.front().id].right.empty())
            leaf_nodes += 1.0;
        for (auto rhs: vector_rules[q_nt1.front().id].right)
            if (rhs.first.size() == 2) {
                bool exist_both_subtree = false;
                for (auto rhs2: vector_rules[q_nt2.front().id].right)
                    if (rhs2.first.size() == 2)
                        if (rhs.first[0].equal_symbol(rhs2.first[0])) {
                            if (rhs.first[1].name.empty())
                                cout << "error here" << endl;
                            n_sons +=1.0;
                            q_nt1.push(rhs.first[1]);
                            q_nt2.push(rhs2.first[1]);
                            exist_both_subtree = true;
                            //if (exist_empty_rule(rules[rhs2.first[1].id].right))
                            if (rules[rhs2.first[1].id].right.size() == 1)
                                exist_both_subtree = false;
                            break;
                        }
                if (!exist_both_subtree) {
                    if (q_nt1.front().id == nt1.id)
                        return false;
                    return true;
                }
            }
        q_nt1.pop();
        q_nt2.pop();
    }
    //cout << nt1.name << " * " << nt2.name << " : " << count/n_nodes << endl;
    /*if (count/n_nodes >= tolerance)
        return true;*/
    return false;
}
std::vector<Rule::Rule> Grammar::Grammar::find_pumping_rule_by_auto_similarity(Symbol::Symbol nt, set<int> &not_search_nts, std::map<int, std::pair<std::vector<std::pair<Symbol::Symbol, int>>, int>> compatible_lists, double p_ratio) {
    //cout << "   Checking nt at " << nt.id << "/" << non_terminals.size() << endl ;
    queue<Symbol::Symbol> queue_symbol;
    queue_symbol.push(nt);
    map<int, vector<tuple<vector<Symbol::Symbol>, vector<Symbol::Symbol>, vector<Symbol::Symbol>, vector<Symbol::Symbol>, set<int>, int>>> nts_pumpings;
    while (!queue_symbol.empty()) {
        for (auto rhs: rules[queue_symbol.front().id].right) {
            if (rhs.first.size() == 2) {
                queue_symbol.push(rhs.first[1]);
            }
        }
        bool flag = exist_empty_rule(rules[queue_symbol.front().id].right);
        /*if (flag)
            cout << "";*/
        if (flag) {

            int frequency = 0;
            int amount = 0;
            vector<Symbol::Symbol> z_path;
            vector<Symbol::Symbol> path_to_accept;
            path_to_accept.push_back(queue_symbol.front());
            while (path_to_accept[path_to_accept.size()-1].id != nt.id) {
                for (auto r: rules) {
                    bool back_nt_found = false;
                    for (auto rhs: r.right) {
                        if (rhs.first.size() == 2) {
                            if (rhs.first[1].id == path_to_accept[path_to_accept.size()-1].id) {
                                path_to_accept.push_back(r.left[0]);
                                back_nt_found = true;
                                break;
                            }
                        }
                    }
                    if (back_nt_found)
                        break;

                }
            }
            //cout  << path_to_accept[0].name << ", ";

            while (path_to_accept.size() > 1) {
                set<int> nts_that_pumps;
                int count_repeatence = 0;
                queue<pair<Symbol::Symbol, int>> queue_symbol_auto_similarity;
                vector<pair<Symbol::Symbol, int>> list_symbol_auto_similarity;
                list_symbol_auto_similarity = compatible_lists[path_to_accept[path_to_accept.size()-1].id].first;
                /*for (auto rhs: rules[nt.id].right)
                    if (rhs.first.size() == 2)
                        queue_symbol_auto_similarity.push(make_pair(rhs.first[1], 1));

                while (!queue_symbol_auto_similarity.empty()) {
                    if (fpta_pumping_compatible_tree2( path_to_accept[path_to_accept.size()-1], queue_symbol_auto_similarity.front().first,1.5, rules, non_terminals)) {//if (check_auto_similarity(path_to_accept, queue_symbol_auto_similarity.front().first)) {
                        //fpta_pumping_compatible_tree2( path_to_accept[path_to_accept.size()-1], queue_symbol_auto_similarity.front().first,1.5, rules, non_terminals);
                        list_symbol_auto_similarity.push_back(queue_symbol_auto_similarity.front());
                        count_repeatence++;
                    }
                    for (auto rhs: rules[queue_symbol_auto_similarity.front().first.id].right)
                        if (rhs.first.size() == 2)
                            queue_symbol_auto_similarity.push(make_pair(rhs.first[1], queue_symbol_auto_similarity.front().second+1));
                    queue_symbol_auto_similarity.pop();

                }*/
                set<int> accept_repetition;
                map<int, int> start_pumping;
                int last_start_pump = 0;
                int next_start_pump = 0;
                /*if (!list_symbol_auto_similarity.empty())
                    cout << endl <<  "      found some auto similarities at " << path_to_accept[0].name;*/
                for (auto s: list_symbol_auto_similarity) {


                    if (start_pumping.find(s.second - next_start_pump) == start_pumping.end()) {
                        start_pumping[s.second - last_start_pump] += 1;
                        next_start_pump = last_start_pump;
                        last_start_pump = s.second;
                    } else {
                        start_pumping[s.second - next_start_pump] += 1;
                    }
                    queue<pair<Symbol::Symbol, int>>s_to_first_accept;
                    s_to_first_accept.push(s);
                    while (!s_to_first_accept.empty()) {
                        for (auto rhs: rules[s_to_first_accept.front().first.id].right) {
                            if (rhs.first.size() == 2) {
                                s_to_first_accept.push(make_pair(rhs.first[1], s_to_first_accept.front().second+1));
                            }
                        }
                        if (exist_empty_rule(rules[s_to_first_accept.front().first.id].right) ) {
                            if (s_to_first_accept.front().second - (int) path_to_accept.size() > 0) {
                                //cout << "Accept in: " << s_to_first_accept.front().second -  s.second << endl;
                                accept_repetition.insert(s_to_first_accept.front().second - (int) path_to_accept.size());
                            }
                            break;
                        }
                        s_to_first_accept.pop();
                    }
                }
                //sort(accept_repetition.begin(),accept_repetition.end());
                if (!list_symbol_auto_similarity.empty() && !accept_repetition.empty()) {

                    map<int, int> count_accept_difference;
                    for (auto e: accept_repetition) {
                        std::pair<std::set<int>::const_iterator,std::set<int>::const_iterator> ret;
                        ret = accept_repetition.equal_range(e);
                        count_accept_difference[*ret.second - *ret.first] += 1;

                    }
                    int max = 0;
                    int pumping_size = 0;
                    if (accept_repetition.size() == 1)
                        pumping_size = *accept_repetition.begin();
                    else {
                        for (auto e: count_accept_difference) {
                            if (e.second > max && e.first > 0) {
                                max = e.second;
                                pumping_size = e.first;
                            }
                        }
                    }
                    if (pumping_size == 0)
                        break;
                    if (path_to_accept.size() >= pumping_size) {
                        //cout << endl <<  "      found some auto similarities with repetition at " << path_to_accept[0].name;
                        max = 0;
                        int v_size = 0;
                        for (auto e: start_pumping) {
                            if (e.second > max) {
                                max = e.second;
                                v_size = e.first;
                            }
                        }
                        if (pumping_size < v_size)
                            break;
                        vector<Symbol::Symbol> u;
                        vector<Symbol::Symbol> v;
                        vector<Symbol::Symbol> w;
                        vector<Symbol::Symbol> x;
                        vector<Symbol::Symbol> z;
                        vector<Symbol::Symbol> null_v;

                        find_v_w_x_z_from_path(path_to_accept, v, w, x, v_size, pumping_size, z_path, z);
                        bool exist_derivation = true;
                        int total_yielded_words = 0;
                        int total_words = 0;
                        double total_yielded_words_p = 0;
                        double total_words_p = 0;
                        int pumpings_use = 0;
                        for (int i = 2; w.size()+z.size()+i*(v.size()+x.size()) < compatible_lists[nt.id].second; i++) {
                        //for (int i = 2; i < 5; i++) {

                            total_words_p += 1.0;
                            if (check_derivation_from_nt_with_v_w_x_z(nt, v, w, x, z, i)) {
                                total_yielded_words_p += 1.0;
                                pumpings_use += i;
                            }
                            for (auto e: list_symbol_auto_similarity) {
                                if (e.second % v_size == 0) {
                                    total_words++;
                                    if (check_reachable_node_from_nt_with_v_w_x_z(e.first, v, w, x, z, i, e.second/v_size)) {
                                        nts_that_pumps.insert(e.first.id);
                                        total_yielded_words++;

                                    }
                                }
                            }
                        }

                        if (total_yielded_words_p/total_words_p > p_ratio) { // e nts_that_pumps vazio
                            cout << endl <<  "          found P-RULE! " << "u "<< queue_symbol.front().name << ", v " << convert_vector_to_string(v) <<  ", w " << convert_vector_to_string(w) << ", x " << convert_vector_to_string(x) << ", z " << convert_vector_to_string(z) << ", p_use: " << pumpings_use;
                            nts_pumpings[queue_symbol.front().id].push_back(make_tuple(v, w, x, z, nts_that_pumps, pumping_size));
                        } else if (total_yielded_words_p > 1.0){
                            cout << endl <<  "          not found P-RULE! " << "u "<< queue_symbol.front().name << ", v " << convert_vector_to_string(v) <<  ", w " << convert_vector_to_string(w) << ", x " << convert_vector_to_string(x) << ", z " << convert_vector_to_string(z) << ", total words: " << total_yielded_words_p <<  ", ratio: " << total_yielded_words_p/total_words_p;
                        }

                    } else
                        break;
                }
                z_path.push_back(path_to_accept[0]);
                path_to_accept.erase(path_to_accept.begin(), path_to_accept.begin()+1);
            }
            //if (!nts_pumpings.empty())
                //cout<< endl << "        in word " << queue_symbol.front().name << ", found pumps" << endl << "  continue words: ";
            /*for (auto p: nts_pumpings[queue_symbol.front().id]) {
                //cout << "v: " << convert_vector_to_string(get<0>(p)) << ", w: " << convert_vector_to_string(get<1>(p)) << ", x: " << convert_vector_to_string(get<2>(p)) <<  " and z: " << convert_vector_to_string(get<3>(p)) << " start in nts: ";
                for (auto nt: get<4>(p)) {
                    //cout << non_terminals[nt].name << " ";

                }
                //cout << endl;
            }*/
        }
        queue_symbol.pop();
    }
    vector<Rule::Rule> return_rules = mount_pumping_rules(nts_pumpings, not_search_nts, nt, compatible_lists[nt.id].second);
    return return_rules;
}
bool Grammar::Grammar::check_auto_similarity(std::vector<Symbol::Symbol> path_to_accept, Symbol::Symbol start) {
    for (int i = path_to_accept.size()-1; i >= 0 ; i--) {
        if (!fpta_pumping_compatible_tree2( path_to_accept[i], start,1.5, rules, non_terminals))
        //if (!compatible_alergia(start, path_to_accept[i], 1, rules))
            return false;
        if (i > 0) {
            bool exist_path = false;
            Symbol::Symbol symbol_path;
            for (auto rhs2: rules[path_to_accept[i].id].right) {
                if (rhs2.first.size() == 2) {
                    if (rhs2.first[1].id == path_to_accept[i-1].id) {
                        symbol_path = rhs2.first[0];
                        break;
                    }
                }
            }
            for (auto rhs: rules[start.id].right) {
                if (rhs.first.size() == 2) {
                    if (rhs.first[0].id == symbol_path.id) {
                        start = rhs.first[1];
                        exist_path = true;
                    }
                }
            }
            if (!exist_path)
                return false;
        }
    }
    return true;
}
void Grammar::Grammar::find_v_w_x_z_from_path(std::vector<Symbol::Symbol> path, std::vector<Symbol::Symbol> &v, std::vector<Symbol::Symbol> &w, std::vector<Symbol::Symbol> &x, int v_size, int pumping_size, std::vector<Symbol::Symbol> z_path, std::vector<Symbol::Symbol> &z) {
    reverse(path.begin(), path.end());
    for (int i = 0 ; i < v_size; i ++) {
        for (auto rhs: rules[path[i].id].right) {
            if (rhs.first.size() == 2) {
                if (rhs.first[1].id == path[i+1].id) {
                    v.push_back(rhs.first[0]);
                    break;
                }
            }
        }
    }
    for (int i = v_size ; i < path.size()-pumping_size+v_size-1; i ++) {
        for (auto rhs: rules[path[i].id].right) {
            if (rhs.first.size() == 2) {
                if (rhs.first[1].id == path[i+1].id) {
                    w.push_back(rhs.first[0]);
                    break;
                }
            }
        }
    }
    for (int i = path.size()-pumping_size+v_size-1 ; i < path.size()-1; i ++) {
        for (auto rhs: rules[path[i].id].right) {
            if (rhs.first.size() == 2) {
                if (rhs.first[1].id == path[i+1].id) {
                    x.push_back(rhs.first[0]);
                    break;
                }
            }
        }
    }
    if (!z_path.empty()) {
        reverse(z_path.begin(), z_path.end());
        for (auto rhs: rules[path[path.size()-1].id].right) {
            if (rhs.first.size() == 2) {
                if (rhs.first[1].id == z_path[0].id) {
                    z.push_back(rhs.first[0]);
                    break;
                }
            }
        }

        for (int i = 0 ; i < z_path.size()-1; i ++) {
            for (auto rhs: rules[z_path[i].id].right) {
                if (rhs.first.size() == 2) {
                    if (rhs.first[1].id == z_path[i+1].id) {
                        z.push_back(rhs.first[0]);
                        break;
                    }
                }
            }
        }
    }

}
bool Grammar::Grammar::check_reachable_node_from_nt_with_v_w_x_z(Symbol::Symbol nt, vector<Symbol::Symbol> &v, vector<Symbol::Symbol> &w, vector<Symbol::Symbol> &x, vector<Symbol::Symbol> &z, int pumping_times, int already_pumped_v) {

    if (already_pumped_v == 0)
        return false;
    //Se não precisar derivar o bombeamento V, apagar fors aninhados abaixo
    for (int i = 0; i < pumping_times- already_pumped_v; i++) {
        for (auto v_e: v) {
            bool exist_rhs = false;
            for (auto rhs: rules[nt.id].right) {
                if (rhs.first.size() == 2) {
                    if (rhs.first[0].id == v_e.id) {
                        nt = rhs.first[1];
                        exist_rhs = true;
                        break;
                    }

                }
            }
            if (!exist_rhs)
                return false;
        }
    }
    for (auto w_e: w) {
        bool exist_rhs = false;
        for (auto rhs: rules[nt.id].right) {

            if (rhs.first.size() == 2) {
                if (rhs.first[0].id == w_e.id) {
                    nt = rhs.first[1];
                    exist_rhs = true;
                    break;
                }
            }
        }
        if (!exist_rhs)
            return false;
    }
    for (int i = 0; i < pumping_times; i++) {
        for (auto x_e: x) {
            bool exist_rhs = false;
            for (auto rhs: rules[nt.id].right) {
                if (rhs.first.size() == 2) {
                    if (rhs.first[0].id == x_e.id) {
                        nt = rhs.first[1];
                        exist_rhs = true;
                        break;
                    }
                }
            }
            if (!exist_rhs)
                return false;
        }
    }
    for (auto z_e: z) {
        bool exist_rhs = false;
        for (auto rhs: rules[nt.id].right) {
            if (rhs.first.size() == 2) {
                if (rhs.first[0].id == z_e.id) {
                    nt = rhs.first[1];
                    exist_rhs = true;
                    break;
                }
            }
        }
        if (!exist_rhs)
            return false;
    }

    return true;
}


bool Grammar::Grammar::check_derivation_from_nt_with_v_w_x_z(Symbol::Symbol nt, vector<Symbol::Symbol> &v, vector<Symbol::Symbol> &w, vector<Symbol::Symbol> &x, vector<Symbol::Symbol> &z, int pumping_times) {

    //Se não precisar derivar o bombeamento V, apagar fors aninhados abaixo
    for (int i = 0; i < pumping_times; i++) {
        for (auto v_e: v) {
            bool exist_rhs = false;
            for (auto rhs: rules[nt.id].right) {
                if (rhs.first.size() == 2) {
                    if (rhs.first[0].id == v_e.id) {
                        nt = rhs.first[1];
                        exist_rhs = true;
                        break;
                    }

                }
            }
            if (!exist_rhs)
                return false;
        }
    }
    for (auto w_e: w) {
        bool exist_rhs = false;
        for (auto rhs: rules[nt.id].right) {

            if (rhs.first.size() == 2) {
                if (rhs.first[0].id == w_e.id) {
                    nt = rhs.first[1];
                    exist_rhs = true;
                    break;
                }
            }
        }
        if (!exist_rhs)
            return false;
    }
    for (int i = 0; i < pumping_times; i++) {
        for (auto x_e: x) {
            bool exist_rhs = false;
            for (auto rhs: rules[nt.id].right) {
                if (rhs.first.size() == 2) {
                    if (rhs.first[0].id == x_e.id) {
                        nt = rhs.first[1];
                        exist_rhs = true;
                        break;
                    }
                }
            }
            if (!exist_rhs)
                return false;
        }
    }
    for (auto z_e: z) {
        bool exist_rhs = false;
        for (auto rhs: rules[nt.id].right) {
            if (rhs.first.size() == 2) {
                if (rhs.first[0].id == z_e.id) {
                    nt = rhs.first[1];
                    exist_rhs = true;
                    break;
                }
            }
        }
        if (!exist_rhs)
            return false;
    }

    if (exist_empty_rule(rules[nt.id].right))
        return true;
    return false;
}
bool Grammar::Grammar::check_v_pumping_use(Symbol::Symbol nt, Symbol::Symbol nt_target, std::vector<Symbol::Symbol> &v, int pumping_times) {
    for (int i = 0; i < pumping_times; i++) {
        for (auto v_e: v) {
            bool exist_rhs = false;
            for (auto rhs: rules[nt.id].right) {
                if (rhs.first.size() == 2) {
                    if (rhs.first[0].id == v_e.id) {
                        nt = rhs.first[1];
                        exist_rhs = true;
                        break;
                    }

                }
            }
            if (!exist_rhs)
                return false;
        }
    }
    if (nt.id == nt_target.id)
        return true;
    return false;
}

std::vector<Rule::Rule> Grammar::Grammar::mount_pumping_rules(std::map<int, std::vector<std::tuple<std::vector<Symbol::Symbol>, std::vector<Symbol::Symbol>, std::vector<Symbol::Symbol>, std::vector<Symbol::Symbol>, std::set<int>, int>>> nts_pumpings, set<int> &not_search_nts, Symbol::Symbol nt, int max_height) {
    vector<tuple<vector<Symbol::Symbol>, vector<Symbol::Symbol>, int>> z_word_counts;
    for (auto e1: nts_pumpings) {
        for (auto e2: nts_pumpings) {
            if (e1.first != e2.first) {
                for (auto z1: e1.second) {
                    for (auto z2: e2.second) {
                        if (equal_word(get<3>(z1), get<3>(z2))) {
                            if (equal_word(get<1>(z1), get<1>(z2))) {
                                bool there_is_z_count = false;
                                for (auto &z_count : z_word_counts) {
                                    if (equal_word(get<3>(z1), get<1>(z_count))) {
                                        if (equal_word(get<1>(z1), get<0>(z_count))) {
                                            get<2>(z_count)++;
                                            there_is_z_count = true;
                                            break;
                                        }
                                    }
                                }
                                if (!there_is_z_count)
                                    z_word_counts.push_back(make_tuple(get<1>(z1), get<3>(z1), 1));
                            }
                        }
                    }
                }
            }
        }
    }
    tuple<vector<Symbol::Symbol>, vector<Symbol::Symbol>, int> max_z_count ;
    get<2>(max_z_count) = 0;
    for (auto z: z_word_counts) {
        if (get<2>(z) > get<2>(max_z_count))
            max_z_count = z;
    }
    vector<Rule::Rule> pumping_rules;
    Symbol::Symbol p_nt = Symbol::Symbol("!" + nt.name, 0, false, false);
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right;
    for (auto e : nts_pumpings) {
        for (auto p: e.second) {
            if (equal_word(get<3>(p), get<1>(max_z_count))) {
                if (equal_word(get<1>(p), get<0>(max_z_count))) {
                    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> rhs;
                    rhs.first.insert(rhs.first.end(), get<0>(p).begin(), get<0>(p).end());
                    rhs.first.push_back(p_nt);
                    rhs.first.insert(rhs.first.end(), get<2>(p).begin(), get<2>(p).end());
                    for (int i = 1; i < get<4>(p).size(); i++) {
                        rhs.second.first += 1;
                    }
                    right.push_back(rhs);
                    rhs.first.clear();
                    rhs.second.first = 0.0;
                    rhs.first.insert(rhs.first.end(), get<1>(p).begin(), get<1>(p).end());
                    for (int i = 1; i < get<4>(p).size(); i++)
                        rhs.second.first += 1;
                    right.push_back(rhs);
                    not_search_nts.merge(get<4>(p));
                }
            }
        }
    }
    group_equal_rhs(right);
    std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> w;
    for (int i = 0; i < right.size(); i ++) {
        bool is_pump = false;
        for (auto  & s : right[i].first) {
            if (s.name.find("!#$") != string::npos) {
                is_pump = true;
                break;
            }
        }
        if (!is_pump) {
            w = right[i];

            right.erase(right.begin()+i);
            break;
        }

    }

    right.push_back(w);
    vector<Symbol::Symbol> lhs;
    lhs.push_back(p_nt);
    pumping_rules.push_back(Rule::Rule(lhs, right));
    lhs.clear();
    right.clear();
    lhs.push_back(nt);
    std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rhs;
    rhs.first.push_back(p_nt);
    rhs.first.insert(rhs.first.end(), get<1>(max_z_count).begin(), get<1>(max_z_count).end());
    rhs.second.first = 0.0;
    right.push_back(rhs);
    pumping_rules.push_back(Rule::Rule(lhs, right));
    if (pumping_rules[0].freq() > 0.0)
        decrease_rules_from_pumping(nt, pumping_rules, max_height, get<1>(max_z_count));
    for (auto  & rhs: pumping_rules[0].right)
        if (rhs.first.empty())
            rhs.first.push_back(Symbol::Symbol("", -1, true, false));
    return pumping_rules;
}
void Grammar::Grammar::super_duper_pumping_inference(double alpha, double p_ratio, double time_limite) {
    std::chrono::duration<double> iterationTime = std::chrono::steady_clock::now() -  std::chrono::steady_clock::now();
    auto startIt = std::chrono::steady_clock::now();
    non_terminals.clear();
    rules.clear();
    n_non_terminals = 0;
    gen_fpta();
    print_grammar();
    //cout <<" Finished FPTA" << endl;
    std::map<int, std::pair<std::vector<std::pair<Symbol::Symbol, int>>, int>> compatible_lists = build_compatible_lists(alpha);
    set<int> not_search_nts;
    vector<Rule::Rule> pumped_rules;
    vector<Rule::Rule> pumping_rules;
    vector<Symbol::Symbol> pumped_nts;
    int count_nt  = 0;
    for (auto nt: non_terminals) {
        if ((count_nt+1) % (n_non_terminals/10) == 0 )
            cout << "   Pumping inference at " << (100.0*(count_nt+1)/(n_non_terminals*1.0)) << "%" << endl;
        if (not_search_nts.find(nt.id) == not_search_nts.end()) {
            vector<Rule::Rule> new_rules = find_pumping_rule_by_auto_similarity(nt, not_search_nts, compatible_lists, p_ratio);

            if (new_rules[0].freq() > 0.0) {
                //pump_and_reduce_w(new_rules[0]);
                new_rules[1].right[0].first[0].id = n_non_terminals;
                n_non_terminals++;
                non_terminals.push_back(new_rules[1].right[0].first[0]);
                pumping_rules.push_back(new_rules[0]);
                pumped_rules.push_back(new_rules[1]);
                pumped_nts.push_back(new_rules[1].left[0]);
            }
        }
        count_nt++;
        if (time_limite < std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startIt).count()/1000) {
            //cout << "   Time Limit Reached" << endl;
            //break;
        }
    }
    for (auto r: pumped_rules) {
        /*if (!verify_nt_in_subtree_if_any_pumped_not_equal(r.left[0], pumped_nts)) {
            rules[r.left[0].id].right = r.right;

        } else {*/
            rules[r.left[0].id].right.insert(rules[r.left[0].id].right.end(), r.right.begin(), r.right.end());
        //}
    }

    vector<Rule::Rule> aux_rules;
    vector<Symbol::Symbol> aux_nts;
    queue<Symbol::Symbol> nt_queue;
    nt_queue.push(non_terminals[0]);

    while (!nt_queue.empty()) {
        aux_rules.push_back(rules[nt_queue.front().id]);
        aux_nts.push_back(nt_queue.front());
        for (auto rhs: rules[nt_queue.front().id].right) {
            if (rhs.first.size() == 2) {
                if (rhs.first[0].name[0] != '!')  {
                    nt_queue.push(rhs.first[1]);
                }
            }
        }
        nt_queue.pop();
    }


    aux_rules.insert(aux_rules.end(), pumping_rules.begin(), pumping_rules.end());
    for (auto pr: pumping_rules)
        aux_nts.push_back( pr.left[0]);
    rules = aux_rules;
    non_terminals = aux_nts;


    for (int i = 0; i < non_terminals.size(); i++)
        non_terminals[i].id = i;
    n_non_terminals = non_terminals.size();

    for (int i = 0; i < rules.size(); i ++) {
        for (int j = 0; j < rules[i].left.size(); j++)
            rules[i].left[j] = find_symbol_in_vector_by_name(rules[i].left[j].name, non_terminals);
        for (int j = 0; j < rules[i].right.size(); j++)
            for (int k = 0; k < rules[i].right[j].first.size(); k++)
                if (!rules[i].right[j].first[k].terminal)
                    rules[i].right[j].first[k] = find_symbol_in_vector_by_name(rules[i].right[j].first[k].name, non_terminals);
    }

    nullify_unreacheble_rules();
    remove_unused_rules();
    for (auto & r: rules)
        remove_unused_rule_zero_righties(r.right);
    generate_nt_for_t();
    normalize_probs();
    iterationTime = std::chrono::steady_clock::now() - startIt;
    cout << "RunTime: " << std::chrono::duration_cast<std::chrono::milliseconds>(iterationTime).count();
    //print_grammar();
}
void Grammar::Grammar::pump_and_reduce_w(Rule::Rule &r) {
    bool found_w = false;
    Rule::Rule aux_r = r;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right;
    vector<vector<Symbol::Symbol>> w_words;
    for (int i = 0; i < r.right.size(); i ++) {
        bool found_p = false;
        for (auto s: r.right[i].first)
            if (s.name[0] == '!')
                found_p = true;
        if (!found_p) {
            if (!found_w) {
                right.push_back(r.right[i]);
                found_w = true;
            }
            else
                w_words.push_back(r.right[i].first);
        } else
            right.push_back(r.right[i]);
    }
    aux_r.right = right;
    Grammar g = Grammar({}, 1, {}, g.pcfg, make_pair(0, 0));
    g.non_terminals.clear();
    g.rules.clear();
    g.terminals = terminals;
    g.rules.push_back(aux_r);
    g.non_terminals.push_back(r.left[0]);
    g.n_non_terminals = 1;


    for (auto t: g.terminals) {
        Symbol::Symbol new_ntt = Symbol::Symbol("#$$"+t.name, g.n_non_terminals, false, false);
        g.non_terminals.push_back(new_ntt);
        g.n_non_terminals++;
        for (auto &r2 : g.rules)
            for (auto & rhs: r2.right) {
                //rhs.second.first = 0.0;
                for (auto & s: rhs.first) {
                    if (s.name.compare(t.name) == 0)
                        s = new_ntt.clone();
                    else if (s.name.compare(new_ntt.name) == 0)
                        s = new_ntt.clone();
                }
            }
        vector<Symbol::Symbol> left;
        left.push_back(new_ntt);
        pair<vector<Symbol::Symbol>,pair<double, double>> right;
        right.first.push_back(t);
        right.second.first = 1.0;
        vector<pair<vector<Symbol::Symbol>,pair<double, double>>> righties;
        righties.push_back(right);
        Rule::Rule r2 = Rule::Rule(left,righties);
        g.rules.push_back(r2);

    }

    for (auto & w: w_words) {
        if (g.calculate_parse_tree_prob_top_down(w) > 0) {
            for (auto & rhs: r.right) {
                if (equal_rhs(rhs.first, w)) {
                    rhs.second.first = 0;
                }
            }
        }
    }
    remove_unused_rule_zero_righties(r.right);


}
std::map<int, std::pair<std::vector<std::pair<Symbol::Symbol, int>>, int>> Grammar::Grammar::build_compatible_lists(double alpha) {
    std::map<int, std::pair<std::vector<std::pair<Symbol::Symbol, int>>, int>> compatible_lists;
    int count_nt = 0;
    for (auto nt: non_terminals) {
        int max_height = 0;
        /*if ((count_nt+1) % (n_non_terminals/10) == 0 )
            cout << "   Build compatible list  at " << (100.0*(count_nt+1)/(n_non_terminals*1.0)) << "%" << endl;*/
        cout << "   Build compatible list  at " << nt.id << "/" << non_terminals.size() << endl ;
        queue<pair<Symbol::Symbol, int>> queue_symbol_auto_similarity;
        vector<pair<Symbol::Symbol, int>> list_symbol_auto_similarity;
        for (auto rhs: rules[nt.id].right)
            if (rhs.first.size() == 2)
                queue_symbol_auto_similarity.push(make_pair(rhs.first[1], 1));
        while (!queue_symbol_auto_similarity.empty()) {
            if (queue_symbol_auto_similarity.front().second > max_height)
                max_height = queue_symbol_auto_similarity.front().second;
            if (fpta_pumping_compatible_tree2( nt, queue_symbol_auto_similarity.front().first,alpha, rules, non_terminals)) {//if (check_auto_similarity(path_to_accept, queue_symbol_auto_similarity.front().first)) {
                //fpta_pumping_compatible_tree2( path_to_accept[path_to_accept.size()-1], queue_symbol_auto_similarity.front().first,1.5, rules, non_terminals);
                list_symbol_auto_similarity.push_back(queue_symbol_auto_similarity.front());
            }
            for (auto rhs: rules[queue_symbol_auto_similarity.front().first.id].right)
                if (rhs.first.size() == 2)
                    queue_symbol_auto_similarity.push(make_pair(rhs.first[1], queue_symbol_auto_similarity.front().second+1));
            queue_symbol_auto_similarity.pop();

        }
        compatible_lists[nt.id] = make_pair(list_symbol_auto_similarity, max_height);
        count_nt++;
    }
    return compatible_lists;
}
void Grammar::Grammar::decrease_rules_from_pumping(Symbol::Symbol nt, std::vector<Rule::Rule> & rs, int max_height, std::vector<Symbol::Symbol> & z) {
    vector<vector<pair<int,int>>> parse_trees;
    vector<pair<int,int>> parse_tree_indexes;
    recursive_build_parse_tree_indexes(rs, max_height, z, parse_tree_indexes, parse_trees);
    for (auto & rhs: rs[0].right)
        rhs.second.first = 0.0;
    for (auto pt: parse_trees) {
        if (check_reacheable_multiple_pumpings(nt, pt, rs))
            rs[1].right[0].second.first +=pt.size()-1;
    }
}
void Grammar::Grammar::recursive_build_parse_tree_indexes(std::vector<Rule::Rule> & rs, int & max_height, std::vector<Symbol::Symbol> & z, std::vector<std::pair<int,int>> parse_tree_indexes, std::vector<std::vector<std::pair<int,int>>> & parse_trees) {
    int tree_height = rs[0].right[rs[0].right.size()-1].first.size() + z.size();
    for (int i = 0; i < parse_tree_indexes.size(); i++)
        tree_height += rs[parse_tree_indexes[i].first].right[parse_tree_indexes[i].second].first.size()-1;

    if (tree_height <= max_height) {


        for (int i = 0; i < rs[0].right.size()-1; i++) {
            std::vector<std::pair<int,int>> aux = parse_tree_indexes;
            aux.push_back(make_pair(0,i));
            recursive_build_parse_tree_indexes(rs, max_height, z, aux, parse_trees);
        }
        parse_tree_indexes.push_back(make_pair(0, rs[0].right.size()-1));
        parse_trees.push_back(parse_tree_indexes);
    }

}
bool Grammar::Grammar::check_reacheable_multiple_pumpings(Symbol::Symbol nt, std::vector<std::pair<int, int>> parse_tree_indexes, vector<Rule::Rule> &rs) {

    vector<Rule::Rule> aux_rules = rules;
    vector<Rule::Rule> aux_p_rules = rs;
    vector<Symbol::Symbol> after_pumps;
    for (auto d : parse_tree_indexes) {
        for (int i = 0; i < rs[d.first].right[d.second].first.size(); i++) {
            Symbol::Symbol s = rs[d.first].right[d.second].first[i];
            bool exist_rhs = false;
            if (s.name.find("!#$") == string::npos) {
                for (int j = 0; j < rules[nt.id].right.size(); j++) {
                    std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rhs = rules[nt.id].right[j];
                    if (rhs.first.size() == 2) {
                        if (rhs.first[0].id == s.id) {

                            if (aux_rules[nt.id].right[j].second.first == 0.0) {
                                //cout << "Rule " << nt.name << "-> " << convert_vector_to_string(aux_rules[nt.id].right[j].first) <<"reached zero" << endl;
                                return false;
                            }
                            aux_rules[nt.id].right[j].second.first -= 1;
                            nt = rhs.first[1];
                            exist_rhs = true;
                            break;
                        }
                    }
                }
                if (!exist_rhs)
                    return false;
            } else {
                aux_p_rules[d.first].right[d.second].second.first += 1;
                vector<Symbol::Symbol> aux;
                aux.insert(aux.end(), rs[d.first].right[d.second].first.begin() + i+1, rs[d.first].right[d.second].first.end());
                reverse(aux.begin(), aux.end());
                after_pumps.insert(after_pumps.end(), aux.begin(), aux.end());
                break;
            }
        }
    }
    reverse(after_pumps.begin(), after_pumps.end());
    for (int i = 0; i < after_pumps.size(); i++) {
            Symbol::Symbol s = after_pumps[i];
            bool exist_rhs = false;
            for (int j = 0; j < rules[nt.id].right.size(); j++) {
                std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> rhs = rules[nt.id].right[j];
                if (rhs.first.size() == 2) {
                    if (rhs.first[0].id == s.id) {
                        aux_rules[nt.id].right[j].second.first -= 1;
                        nt = rhs.first[1];
                        exist_rhs = true;
                        break;
                    }
                }
            }
            if (!exist_rhs)
                return false;
    }
    rules = aux_rules;
    aux_p_rules[0].right[aux_p_rules[0].right.size()-1].second.first += 1.0;
    rs = aux_p_rules;
    return true;
}
void Grammar::Grammar::nullify_unreacheble_rules() {
    queue<Symbol::Symbol> queue_nt;
    set<int> used_rules;
    queue_nt.push(non_terminals[0]);
    used_rules.insert(non_terminals[0].id);
    while (!queue_nt.empty()) {
        //used_rules.insert(queue_nt.front().id);
        for (auto rhs: rules[queue_nt.front().id].right) {
            if (rhs.first.size() >= 1) {
                if (!rhs.first[0].name.empty()) {
                    used_rules.insert(queue_nt.front().id);
                    if (rhs.first[0].name[0] != '!') {
                        if (rhs.second.first > 0.0) {
                            used_rules.insert(rhs.first[1].id);
                            queue_nt.push(rhs.first[1]);
                        }
                    } else {
                        used_rules.insert(rhs.first[0].id);
                    }

                }

            }
        }
        queue_nt.pop();
    }
    for (auto &  r: rules) {
        if (used_rules.find(r.left[0].id) == used_rules.end()) {
            for (auto & rhs : r.right)
                rhs.second.first = 0.0;
        }
    }
}
std::vector<std::vector<Symbol::Symbol>> Grammar::Grammar::generate_max_size_words_from_rules(int max_t) {
    std::vector<vector<Symbol::Symbol>> terminal_words;
    stack<vector<Symbol::Symbol>> stack_derivations;
    stack_derivations.push(rules[0].left);
    int repeated_words_n = 0;
    while (!stack_derivations.empty()) {
        vector<Symbol::Symbol> d_aux = stack_derivations.top();
        stack_derivations.pop();
        bool terminal_word = true;
        for (int i = 0; i < d_aux.size(); i++) {
            if (!d_aux[i].terminal) {
                terminal_word = false;
                for (auto rhs: rules[d_aux[i].id].right) {
                    std::vector<Symbol::Symbol> aux_derivation;
                    aux_derivation.insert(aux_derivation.end(), d_aux.begin(), d_aux.begin()+i);
                    if (!rhs.first[0].name.empty())
                        aux_derivation.insert(aux_derivation.end(), rhs.first.begin(), rhs.first.end());
                    aux_derivation.insert(aux_derivation.end(), d_aux.begin()+i+1, d_aux.end());
                    int n_t = 0;
                    for (auto s: aux_derivation)
                        if (s.terminal)
                            n_t++;
                    if (n_t <= max_t)
                        stack_derivations.push(aux_derivation);

                }
                break;
            }
        }
        if (terminal_word) {
            bool exist_word = false;
            for (auto w: terminal_words)
                if (equal_word(w, d_aux)) {
                    exist_word = true;
                    repeated_words_n++;
                    break;
                }
            if (!exist_word) {
                terminal_words.push_back(d_aux);
                cout << "stack_size: " << stack_derivations.size() << ", words: " << terminal_words.size() << ", repeateds: " << repeated_words_n << endl;
            }
        }

    }
    return terminal_words;
}
void Grammar::Grammar::generate_nt_for_t() {
    for (auto t: terminals) {
        Symbol::Symbol new_ntt = Symbol::Symbol("#$$"+t.name, n_non_terminals, false, false);
        non_terminals.push_back(new_ntt);
        n_non_terminals++;
        for (auto & r: rules)
            for (auto & rhs: r.right) {
                for (auto & s: rhs.first) {
                    if (s.name.compare(t.name) == 0)
                        s = new_ntt.clone();
                    else if (s.name.compare(new_ntt.name) == 0)
                        s = new_ntt.clone();
                }
            }
        vector<Symbol::Symbol> left;
        left.push_back(new_ntt);
        pair<vector<Symbol::Symbol>,pair<double, double>> right;
        right.first.push_back(t);
        right.second.first = 1.0;
        vector<pair<vector<Symbol::Symbol>,pair<double, double>>> righties;
        righties.push_back(right);
        Rule::Rule r = Rule::Rule(left,righties);
        rules.push_back(r);
    }
}


double Grammar::Grammar::find_word_probabilities(std::vector<Symbol::Symbol> word) {
    map<vector<pair<int,int>>, map<int, double>> probs;
    double prob = 1.0;
    double sum_prob = 0.0;
    if (word.empty()) {
        for (auto & rhs : rules[0].right)
            for (auto & s: rhs.first)
                if (s.name.empty()) {
                    break;
                }
        return rules[0].right[rules[0].right.size()-1].second.first;
    }
    int right = 0;
    stack<pair<Symbol::Symbol, int>> symbol_rights_stack;
    symbol_rights_stack.push(make_pair(rules[0].left[0], 0));
    vector<Symbol::Symbol> yielded_word;
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>> parse_tree;
    vector<pair<int,int>> parse_tree_indexes;
    while (!symbol_rights_stack.empty() && symbol_rights_stack.top().second < rules[symbol_rights_stack.top().first.id].right.size() || !symbol_rights_stack.empty() && symbol_rights_stack.top().second < rules[symbol_rights_stack.top().first.id].right.size() && symbol_rights_stack.top().first.terminal) {

        if (!symbol_rights_stack.top().first.terminal) {
            Symbol::Symbol to_yield = symbol_rights_stack.top().first;
            right = symbol_rights_stack.top().second;
            symbol_rights_stack.pop();

            for (int i = rules[to_yield.id].right[right].first.size()-1; i >=0; i--)
                symbol_rights_stack.push(make_pair(rules[to_yield.id].right[right].first[i], 0));

            parse_tree_indexes.push_back(make_pair(to_yield.id, right));

        } else {
            yielded_word.push_back(symbol_rights_stack.top().first);

            bool flag_yield_ok = true;
            vector<Symbol::Symbol> yielded_no_emptystr_word = remove_empty_substring(yielded_word);

            for (int j = 0; j < yielded_no_emptystr_word.size(); j++) {
                if (!yielded_no_emptystr_word[j].equal_symbol(word[j]) || yielded_no_emptystr_word.size() > word.size()) {
                    flag_yield_ok = false;
                    while (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1 == parse_tree_indexes[parse_tree_indexes.size()-1].second) {
                        if (parse_tree_indexes.size() == 1 && parse_tree_indexes[parse_tree_indexes.size()-1].second == rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1)
                            break;
                        if (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first[0].terminal) {
                            if (!symbol_rights_stack.top().first.equal_symbol(yielded_word[yielded_word.size()-1]))
                                symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size()-1], 0));
                            yielded_word.pop_back();
                        }
                        for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                            symbol_rights_stack.pop();
                        symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0],parse_tree_indexes[parse_tree_indexes.size()-1].second));
                        //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                        parse_tree_indexes.pop_back();
                    }
                    if (!yielded_word.empty()) {
                        if (yielded_word[yielded_word.size()-1].name.empty()) {
                            symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size() - 1], 0));
                            yielded_word.pop_back();
                        }
                    }
                    if (!parse_tree_indexes.empty()) {
                        for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                            symbol_rights_stack.pop();
                        symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0], parse_tree_indexes[parse_tree_indexes.size()-1].second+1));
                        //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                        parse_tree_indexes.pop_back();
                    }
                    break;
                }
            }
            if (flag_yield_ok) {
                symbol_rights_stack.pop();
            }

        }
        if (symbol_rights_stack.empty()) {
            if (equal_word(remove_empty_substring(yielded_word), word)) {
                for (auto p: parse_tree_indexes)
                    prob *= rules[p.first].right[p.second].second.first;
                sum_prob += prob;
                prob = 1.0;
            }

                while (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1 == parse_tree_indexes[parse_tree_indexes.size()-1].second) {
                    if (parse_tree_indexes.size() == 1 && parse_tree_indexes[parse_tree_indexes.size()-1].second == rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1)
                        break;
                    if (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first[0].terminal) {
                        symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size()-1], 0));
                        yielded_word.pop_back();
                    }
                    for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                        symbol_rights_stack.pop();
                    symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0],parse_tree_indexes[parse_tree_indexes.size()-1].second));
                    //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                    parse_tree_indexes.pop_back();
                }
                if (!yielded_word.empty()) {
                    if (yielded_word[yielded_word.size() - 1].name.empty()) {
                        symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size() - 1], 0));
                        yielded_word.pop_back();
                    }
                }
                if (!parse_tree_indexes.empty()) {
                    for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                        symbol_rights_stack.pop();
                    symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0], parse_tree_indexes[parse_tree_indexes.size()-1].second+1));
                    //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                    parse_tree_indexes.pop_back();
                }

        }
        if (symbol_rights_stack.size() == 1 && symbol_rights_stack.top().first.name.empty()) {

            yielded_word.push_back(symbol_rights_stack.top().first);
            symbol_rights_stack.pop();
            if (equal_word(remove_empty_substring(yielded_word), word)) {
                for (auto p : parse_tree_indexes)
                    prob *= rules[p.first].right[p.second].second.first;
                sum_prob += prob;
                prob = 1.0;
            }
            while (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1 == parse_tree_indexes[parse_tree_indexes.size()-1].second) {
                if (parse_tree_indexes.size() == 1 && parse_tree_indexes[parse_tree_indexes.size()-1].second == rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right.size()-1)
                    break;
                if (rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first[0].terminal) {
                    symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size()-1], 0));
                    yielded_word.pop_back();
                }
                for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                    symbol_rights_stack.pop();
                symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0],parse_tree_indexes[parse_tree_indexes.size()-1].second));
                //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                parse_tree_indexes.pop_back();
            }
            if (!yielded_word.empty()) {
                if (yielded_word[yielded_word.size() - 1].name.empty()) {
                    symbol_rights_stack.push(make_pair(yielded_word[yielded_word.size() - 1], 0));
                    yielded_word.pop_back();
                }
            }
            if (!parse_tree_indexes.empty()) {
                for (int k = 0; k < rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].first.size(); k++)
                    symbol_rights_stack.pop();
                symbol_rights_stack.push(make_pair(rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].left[0], parse_tree_indexes[parse_tree_indexes.size()-1].second+1));
                //rules[parse_tree_indexes[parse_tree_indexes.size()-1].first].right[parse_tree_indexes[parse_tree_indexes.size()-1].second].second.first -= 1.0;
                parse_tree_indexes.pop_back();
            }

            //probs[parse_tree_indexes][-1] = 1.0;
            //break;
        }
        while (symbol_rights_stack.top().first.name.empty()) {
            yielded_word.push_back(symbol_rights_stack.top().first);
            symbol_rights_stack.pop();
        }

    }
    return sum_prob;
}
double Grammar::Grammar::find_word_probabilities_from_pcfg_inside_table(std::vector<Symbol::Symbol> word) {
    double ***iTable  = cyk_prob_vec(word);
    double pX = iTable[0][0][word.size()-1] ;
    free_inside_table(iTable, word.size());
    return pX;
}
