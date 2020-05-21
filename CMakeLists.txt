cmake_minimum_required(VERSION 3.14)
project(grammatical-inferences)

set(CMAKE_CXX_STANDARD 14)

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
    Symbol.h)

target_link_libraries(grammatical-inferences stdc++fs)
