//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_SYMBOL_H
#define GRAMMARINDCTION_SYMBOL_H
#include <string>
#include <vector>


class Symbol{
public:
    std::string name;
    unsigned int id{};
    bool terminal{};
    bool context{};
public:
    Symbol(std::string name, unsigned int id, bool terminal);
    Symbol(std::string name, unsigned int id, bool terminal, bool context);

    Symbol();

    [[nodiscard]] Symbol clone() const;
    void printSymbol() const;
    [[nodiscard]] bool equalSymbol(const Symbol& s) const;

};

#endif //GRAMMARINDCTION_SYMBOL_H
