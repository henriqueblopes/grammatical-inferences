//
// Created by henrique on 02/04/19.
//

#include <iostream>
#include <cmath>
#include "Rule.h"

Rule::Rule(const std::vector<Symbol> &left, const std::vector<std::pair<std::vector<Symbol>, std::pair<double, double>>> &right) : left(left),
                                                                                                                right(right) {}

void Rule::printRule() {
    std::vector<std::pair<std::vector<Symbol>, std::pair<double, double>>>::iterator itVector;
    std::vector<Symbol>::iterator itRule;

    //std::cout << "P_D_A: " << probDirichletTheta << " ";
    for (itRule = left.begin(); itRule != left.end(); itRule++)
        std::cout << (*itRule).name << " ";
    std::cout << " <- ";

    for (itVector = right.begin(); itVector != right.end(); itVector++) {
        if ((*itVector).second.first > 0.0001) {
            std::cout << "P: " << (*itVector).second.first << " ";//  << " a: " << (*itVector).second.second << " ";
            for(itRule = (*itVector).first.begin(); itRule != (*itVector).first.end(); itRule++)
                std::cout <<  (*itRule).name << " ";
            if (itVector != right.end()-1)
                std::cout << " / ";
        }

    }
    std::cout << std::endl;

}

std::string Rule::ruleToStr() {
    std::string ruleStr = R"()";
    std::vector<std::pair<std::vector<Symbol>, std::pair<double, double>>>::iterator itVector;
    std::vector<Symbol>::iterator itRule;

    for (itRule = left.begin(); itRule != left.end(); itRule++)
        ruleStr += (*itRule).name;
    ruleStr += " <- ";

    for (itVector = right.begin(); itVector != right.end(); itVector++) {
        for(itRule = (*itVector).first.begin(); itRule != (*itVector).first.end(); itRule++) {
            if ((*itRule).terminal == true)
                ruleStr += "'"+ (*itRule).name + "'"+" ";
            else
                ruleStr += (*itRule).name + " ";
        }

        if (itVector != right.end()-1)
            ruleStr += " / ";
    }
    return ruleStr;

}

std::string Rule::ruleToStrLALR() {
    std::string ruleStr = "";

    std::vector<std::pair<std::vector<Symbol>, std::pair<double, double>>>::iterator itVector;
    std::vector<Symbol>::iterator itRule;

    for (itRule = left.begin(); itRule != left.end(); itRule++)
        ruleStr += (*itRule).name;
    ruleStr += " <- ";

    for (itVector = right.begin(); itVector != right.end(); itVector++) {
        for(itRule = (*itVector).first.begin(); itRule != (*itVector).first.end(); itRule++) {
            if ((*itRule).terminal == true)
                ruleStr += "'"+ (*itRule).name + "'"+" ";
            else
                ruleStr += (*itRule).name + " ";
        }
        std::string s = std::to_string((*itVector).second.first);
        ruleStr += s;
        ruleStr += " ";

        if (itVector != right.end()-1)
            ruleStr += "| ";
        else {
            ruleStr.pop_back();
            ruleStr += ";\n";
        }
    }
    return ruleStr;

}



Rule Rule::clone() {
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>> ps;
    for (std::pair<std::vector<Symbol>,std::pair<double, double>> p : right) {
        std::pair<std::vector<Symbol>,std::pair<double, double>> pc;
        pc.second = p.second;
        for (Symbol s: p.first)
            pc.first.push_back(s.clone());
        ps.push_back(pc);
    }
    std::vector<Symbol> l;
    for (Symbol s : left)
        l.push_back(s.clone());
    return Rule(l,ps);
}

void Rule::generatePiorDirichlet(double alfa, int thetaDist) {
    //alfa = 10E-4 and theta uniform suggested
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
    double theta = 1.0/right.size();
    for (itRight = right.begin(); itRight != right.end(); itRight++) {
        (*itRight).second.first = theta;
        (*itRight).second.second = alfa;
    }
}

double Rule::cConstant() {
    double numerator = 1.0;
    double denominator = 0.0;
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = right.begin(); itRight != right.end(); itRight++) {
        numerator = numerator * tgamma((*itRight).second.second);
        denominator = denominator + (*itRight).second.second;
    }
    return numerator/tgamma(denominator);
}

void Rule::updateProbDirichletTheta() {
    std::vector<std::pair<std::vector<Symbol>,std::pair<double, double>>>::iterator itRight;
    double prodThetaAlfa = 1.0;
    for (itRight = right.begin(); itRight != right.end(); itRight++) {
        prodThetaAlfa = prodThetaAlfa * pow((*itRight).second.first, (*itRight).second.second-1);
    }
    //TO DO: talvez colocar cConstant como atributo
    probDirichletTheta = prodThetaAlfa/cConstant();
}

std::pair<std::vector<Symbol>, std::pair<double, double>> Rule::getRightSidebyId(int id1stNonContext, int id2ndNonContext, bool terminal, int nNonterminals) {
    if (terminal)
        return right[nNonterminals*nNonterminals + id1stNonContext];
    else
        return right[nNonterminals*id1stNonContext + id2ndNonContext];
}



