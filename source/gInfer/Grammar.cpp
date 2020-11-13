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
    std::cout << "Terminals: " << std::endl;
    std::vector<Symbol::Symbol>::iterator it;
    for(it = terminals.begin(); it != terminals.end(); it++) {
        std::cout << "\t";
        std::cout <<  (*it).name << std::endl;
    }

    std::cout << std::endl;

    std::cout << "NonTerminals: " << std::endl;
    for(it = non_terminals.begin(); it != non_terminals.end(); it++) {
        std::cout << "\t";
        std::cout <<  (*it).name <<std::endl;
    }

    std::cout << std::endl;

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
                        //TO DO: Colocar o FOR abaixo dentro do anterior Unindo os pares de não terminais com terminais
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


/*Grammar::Grammar(Grammar const &grammar) {

}*/



