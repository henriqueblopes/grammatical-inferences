//
// Created by henrique on 02/04/19.
//
#include <iostream>
#include <utility>
#include "Symbol.h"


[[maybe_unused]] void Symbol::Symbol::print_symbol() const {
    std::cout << "nome: " << name << " id: " << id;

}

Symbol::Symbol::Symbol(std::string name, size_t id, bool terminal) : name(std::move(name)), id(id), terminal(terminal) {context = false;}

Symbol::Symbol Symbol::Symbol::clone() const {
    return Symbol (name,id,terminal,context);
}

Symbol::Symbol::Symbol(std::string name, size_t id, bool terminal, bool context) : name(std::move(name)), id(id), terminal(terminal), context(context) {

}

bool Symbol::Symbol::equal_symbol(const Symbol& s) const {
    return name == s.name;
}


Symbol::Symbol::Symbol() = default;
