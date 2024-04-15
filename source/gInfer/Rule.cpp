//
// Created by henrique on 02/04/19.
//

#include <iostream>
#include <cmath>
#include <utility>
#include "Rule.h"

Rule::Rule::Rule(std::vector<Symbol::Symbol> left, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right) : left(std::move(left)),
                                                                                                                right(std::move(right)) {}

void Rule::Rule::print_rule() {
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>::iterator itVector;
    std::vector<Symbol::Symbol>::iterator itRule;

    //std::cout << "P_D_A: " << prob_dirichlet_theta << " ";
    for (itRule = left.begin(); itRule != left.end(); itRule++)
        std::cout << (*itRule).id << " - " << (*itRule).name << " ";
    std::cout << " <- ";
    double totalP = 0.0;
    for (itVector = right.begin(); itVector != right.end(); itVector++) {
        totalP += (*itVector).second.first;
        if ((*itVector).second.first > 0.0001) {
            std::cout << "P: " << (*itVector).second.first << " - ";//  << " a: " << (*itVector).second.second << " ";
            for(itRule = (*itVector).first.begin(); itRule != (*itVector).first.end(); itRule++)
                std::cout <<  (*itRule).name << " ";
            if (itVector != right.end()-1)
                std::cout << " / ";
        }

    }
    std::cout << " / " << totalP;
    std::cout << std::endl;

}

std::string Rule::Rule::rule_to_str() {
    std::string ruleStr;
    ruleStr = R"()";
    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>::iterator itVector;
    std::vector<Symbol::Symbol>::iterator itRule;

    for (itRule = left.begin(); itRule != left.end(); itRule++)
        ruleStr += (*itRule).name;
    ruleStr += " <- ";

    for (itVector = right.begin(); itVector != right.end(); itVector++) {
        for(itRule = (*itVector).first.begin(); itRule != (*itVector).first.end(); itRule++) {
            if ((*itRule).terminal)
                ruleStr += "'"+ (*itRule).name + "'"+" ";
            else
                ruleStr += (*itRule).name + " ";
        }

        if (itVector != right.end()-1)
            ruleStr += " / ";
    }
    return ruleStr;

}

std::string Rule::Rule::rule_to_str_lalr() {
    std::string ruleStr;

    std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>>::iterator itVector;
    std::vector<Symbol::Symbol>::iterator itRule;

    for (itRule = left.begin(); itRule != left.end(); itRule++)
        ruleStr += (*itRule).name;
    ruleStr += " <- ";

    for (itVector = right.begin(); itVector != right.end(); itVector++) {
        for(itRule = (*itVector).first.begin(); itRule != (*itVector).first.end(); itRule++) {
            if ((*itRule).terminal)
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



Rule::Rule Rule::Rule::clone() {
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> ps;
    for (const std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>& p : right) {
        std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> pc;
        pc.second = p.second;
        for (const Symbol::Symbol& s: p.first)
            pc.first.push_back(s.clone());
        ps.push_back(pc);
    }
    std::vector<Symbol::Symbol> l;
    for (const Symbol::Symbol& s : left)
        l.push_back(s.clone());
    return Rule(l,ps);
}

void Rule::Rule::generate_pior_dirichlet(double alfa) {
    //alfa = 10E-4 and theta uniform suggested
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    double theta = 1.0/right.size();
    for (itRight = right.begin(); itRight != right.end(); itRight++) {
        (*itRight).second.first = theta;
        (*itRight).second.second = alfa;
    }
}

double Rule::Rule::c_constant() {
    double numerator = 1.0;
    double denominator = 0.0;
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    for (itRight = right.begin(); itRight != right.end(); itRight++) {
        numerator = numerator * tgamma((*itRight).second.second);
        denominator = denominator + (*itRight).second.second;
    }
    return numerator/tgamma(denominator);
}

void Rule::Rule::update_prob_dirichlet_theta() {
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>>::iterator itRight;
    double prodThetaAlfa = 1.0;
    for (itRight = right.begin(); itRight != right.end(); itRight++) {
        prodThetaAlfa = prodThetaAlfa * pow((*itRight).second.first, (*itRight).second.second-1);
    }
    //TO DO: talvez colocar c_constant como atributo
    prob_dirichlet_theta = prodThetaAlfa / c_constant();
}

std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>> Rule::Rule::get_right_side_by_id(size_t id_1st_non_context, size_t id_2nd_non_context, bool terminal, size_t n_non_terminals) {
    if (terminal)
        return right[n_non_terminals * n_non_terminals + id_1st_non_context];
    else
        return right[n_non_terminals * id_1st_non_context + id_2nd_non_context];
}

double Rule::Rule::freq() {
    double freq = 0.0;
    for (const auto& a: right) {
        freq += a.second.first;
    }
    return freq;
}



