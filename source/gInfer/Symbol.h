//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_SYMBOL_H
#define GRAMMARINDCTION_SYMBOL_H
#include <string>
#include <vector>

namespace Symbol {
class Symbol{
public:
    std::string name;
    size_t id{};
    bool terminal{};
    bool context{};
public:
    Symbol(std::string name, unsigned int id, bool terminal);
    Symbol(std::string name, unsigned int id, bool terminal, bool context);

    Symbol();

    [[nodiscard]] Symbol clone() const;
    [[maybe_unused]] void print_symbol() const;
    [[nodiscard]] bool equal_symbol(const Symbol& s) const;

};
}
#endif //GRAMMARINDCTION_SYMBOL_H
