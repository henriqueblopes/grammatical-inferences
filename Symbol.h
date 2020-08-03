//
// Created by henrique on 02/04/19.
//

#ifndef GRAMMARINDCTION_SYMBOL_H
#define GRAMMARINDCTION_SYMBOL_H
#include <string>
#include <vector>


class Symbol{
public:
    Symbol(const std::vector<Symbol> vector, std::pair<int, int> pair, int i, std::vector<std::vector<Symbol>> vector1,
           int i1);

    std::string name;
    unsigned int id;
    bool terminal;
    bool context;
public:
    Symbol(const std::string &name, unsigned int id, bool terminal);
    Symbol(const std::string &name, unsigned int id, bool terminal, bool context);

    Symbol();

    Symbol clone();
    void printSymbol();
    bool equalSymbol(Symbol s);
    static int xplusy (int x, int y);

protected:
    static inline void free(Symbol& s);
};

#endif //GRAMMARINDCTION_SYMBOL_H
