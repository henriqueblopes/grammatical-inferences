//
// Created by henrique on 19/04/25.
//


#ifndef GRAMMATICAL_INFERENCES_RULEH_H
#define GRAMMATICAL_INFERENCES_RULEH_H
#include "Symbol.h"
#include <unordered_map>
#include <vector>

struct VectorHasher {
    int operator()(const std::vector<Symbol::Symbol> &V) const {
        int hash = V.size();
        for(auto &i : V) {
            if (!i.terminal)
                hash ^= i.id + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            else
                hash ^= -(i.id + 0x9e3779b9 + (hash << 6) + (hash >> 2));
        }
        return hash;
    }
};

namespace RuleH {
    class RuleH {
    public:
        std::vector<Symbol::Symbol> left;
        std::unordered_map<std::vector<Symbol::Symbol>,double, VectorHasher> right;

        RuleH(const std::vector<Symbol::Symbol> &left, const std::unordered_map<std::vector<Symbol::Symbol>, double, VectorHasher> &right);
    };
}
#endif//GRAMMATICAL_INFERENCES_RULEH_H
