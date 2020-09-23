//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_RULE_H
#define GRAMMARINDCTION_RULE_H
#include <vector>
#include "Symbol.h"


class Rule {
public:
    std::vector<Symbol> left;
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>> right;
    double probDirichletTheta{};
    int index1stNonContext{};
    std::vector<Symbol> leftContext;
    std::vector<Symbol> rightContext;
    std::vector<std::pair<double, int>> ruleFrequence;
    int index{};


public:
    Rule(std::vector<Symbol> left, std::vector<std::pair<std::vector<Symbol>, std::pair<double, double>>> right);
    void printRule();
    std::string ruleToStr();
    std::string ruleToStrLALR();
    Rule clone();
    void updateProbDirichletTheta();
    void generatePiorDirichlet(double alfa);
    std::pair<std::vector<Symbol>,std::pair<double, double>> getRightSidebyId(int id1stNonContext, int id2ndNonContext, bool terminal, int nNonterminals);
    double freq();
private:
    double cConstant();
};


#endif //GRAMMARINDCTION_RULE_H
