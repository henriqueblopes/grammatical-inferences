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
    double probDirichletTheta;


public:
    Rule(const std::vector<Symbol>  &left, const std::vector<std::pair<std::vector<Symbol>, std::pair<double, double>>> &right);
    void printRule();
    std::string ruleToStr();
    std::string ruleToStrLALR();
    Rule clone();
    void updateProbDirichletTheta();
    void generatePiorDirichlet(double alfa,int thetaDist);
private:
    double cConstant();
};


#endif //GRAMMARINDCTION_RULE_H
