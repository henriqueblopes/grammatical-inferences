//
// Created by henrique on 02/04/19.
//
#include <iostream>
#include <utility>
#include "Symbol.h"


void Symbol::printSymbol() const {
    std::cout << "nome: " << name << " id: " << id;

}

Symbol::Symbol(std::string name, unsigned int id, bool terminal) : name(std::move(name)), id(id), terminal(terminal) {context = false;}

Symbol Symbol::clone() const {
    return Symbol (name,id,terminal,context);
}

Symbol::Symbol(std::string name, unsigned int id, bool terminal, bool context) : name(std::move(name)), id(id), terminal(terminal), context(context) {

}

bool Symbol::equalSymbol(const Symbol& s) const {
    return name == s.name;
}


Symbol::Symbol() = default;
