cmake_minimum_required(VERSION 3.14)
project(grammatical-inferences)

set(CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")





include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(grammatical-inferences
    catch.hpp
    Grammar.cpp
    Grammar.h
    InputWords.cpp
    InputWords.h
    main.cpp
    peglib.h
    Rule.cpp
    Rule.h
    Symbol.cpp
    Symbol.h pybinding.cpp)


find_package(pybind11 REQUIRED)

pybind11_add_module(klpcsg
        Grammar.cpp
        Grammar.h
        Rule.cpp
        Rule.h
        Symbol.cpp
        Symbol.h
        pybinding.cpp)

target_link_libraries(grammatical-inferences stdc++fs ${PYTHON_LIBRARIES} pybind11::embed)
