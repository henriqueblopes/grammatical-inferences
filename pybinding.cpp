//
// Created by henrique on 29/07/2020.
//
#include "Grammar.h"
#include <vector>
#include "Symbol.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


using namespace std;
namespace py = pybind11;
int square(int x) {
    return x * x;
}



PYBIND11_MODULE(klpcsg, m) {

    py::class_<Symbol>(m, "Symbol")
            .def(py::init<const std::string &, unsigned int, bool, bool &>())
            .def("printSymbol", &Symbol::printSymbol)
            .def_readonly("name", &Symbol::name);

    py::class_<Rule>(m, "Rule")
            .def(py::init<const std::vector<Symbol> &, const std::vector<std::pair<std::vector<Symbol>, std::pair<double, double>>> &>())
            .def_readwrite("left", &Rule::left)
            .def_readwrite("right", &Rule::right);

    py::class_<Grammar>(m, "Grammar")
            .def(py::init<const std::vector<Symbol> &, std::pair<int, int>, int, std::vector<std::vector<Symbol>>, int>())
            .def_readwrite("rules", &Grammar::rules)
            .def("train", &Grammar::train)
            .def("grammarToStr", &Grammar::grammarToStr)
            .def("perplexity", &Grammar::perplexityKL);
}