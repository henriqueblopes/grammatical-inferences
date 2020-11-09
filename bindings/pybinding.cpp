//
// Created by henrique on 29/07/2020.
//
#include <gInfer/Grammar.h>
#include <vector>
#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>



using namespace std;
namespace py = pybind11;

//teste
PYBIND11_MODULE(pyGInfer, m) {

    py::class_<Symbol::Symbol>(m, "Symbol")
            .def(py::init<const std::string &, unsigned int, bool, bool &>())
            .def("print_symbol", &Symbol::Symbol::print_symbol)
            .def_readonly("name", &Symbol::Symbol::name);

    py::class_<Grammar::Grammar>(m, "Grammar")
            .def(py::init<const std::vector<Symbol::Symbol> &, int, std::vector<std::vector<Symbol::Symbol>>, Grammar::Grammar::grammar_type, std::pair<int,int>>())
            .def_readwrite("rules", &Grammar::Grammar::rules)
            .def("train", &Grammar::Grammar::train)
            .def("grammar_to_str", &Grammar::Grammar::grammar_to_str)
            .def("perplexity", &Grammar::Grammar::perplexity_kl);
}