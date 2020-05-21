//
// Created by henrique on 02/04/19.
//
#include <iostream>
#include "Symbol.h"
#include "catch.hpp"

void Symbol::free(Symbol &s) {
    //free(s.name);
    //free(s.id);
}

void Symbol::printSymbol() {
    std::cout << "nome: " << name << " id: " << id;

}

Symbol::Symbol(const std::string &name, unsigned int id, bool terminal) : name(name), id(id), terminal(terminal) {}

Symbol Symbol::clone() {
    return Symbol (name,id,terminal,context);
}

Symbol::Symbol(const std::string &name, unsigned int id, bool terminal, bool context) : name(name), id(id), terminal(terminal), context(context) {

}

bool Symbol::equalSymbol(Symbol s) {
    return(!name.compare(s.name));
}

int Symbol::xplusy(int x, int y) {
    return x+y;
}

Symbol::Symbol() {}

TEST_CASE("Test Catch2 for Symbol") {
    REQUIRE (Symbol::xplusy(5,6) == 10);
}