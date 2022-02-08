//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_RULE_H
#define GRAMMARINDCTION_RULE_H
#include <vector>
#include "Symbol.h"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/access.hpp>

namespace Rule{
class Rule {
public:
    std::vector<Symbol::Symbol> left;
    std::vector<std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>>> right;
    double prob_dirichlet_theta{};
    size_t index_1st_non_context{};
    std::vector<Symbol::Symbol> left_context;
    std::vector<Symbol::Symbol> right_context;
    std::vector<std::pair<double, int>> rule_frequence;
    int index{};


    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive & ar, const unsigned int version)
            {
                ar & left;
                ar & right;
                ar & prob_dirichlet_theta;
                ar & index_1st_non_context;
                ar & left_context;
                ar & right_context;
                ar & rule_frequence;
                ar & index;
            }


public:
    Rule(){};
    Rule(std::vector<Symbol::Symbol> left, std::vector<std::pair<std::vector<Symbol::Symbol>, std::pair<double, double>>> right);
    void print_rule();
    std::string rule_to_str();
    std::string rule_to_str_lalr();
    Rule clone();
    void update_prob_dirichlet_theta();
    void generate_pior_dirichlet(double alfa);
    std::pair<std::vector<Symbol::Symbol>,std::pair<double, double>> get_right_side_by_id(size_t id_1st_non_context, size_t id_2nd_non_context, bool terminal, size_t n_non_terminals);
    double freq();
private:
    double c_constant();
};
}

#endif //GRAMMARINDCTION_RULE_H
