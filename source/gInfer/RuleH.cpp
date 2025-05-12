//
// Created by henrique on 19/04/25.
//


#include "RuleH.h"


RuleH::RuleH::RuleH(const std::vector<Symbol::Symbol> &left, const std::unordered_map<std::vector<Symbol::Symbol>, double, VectorHasher> &right) : left(left), right(right) {}
