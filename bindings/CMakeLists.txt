find_package(pybind11 REQUIRED)


pybind11_add_module(pyGInfer
        ${CMAKE_SOURCE_DIR}/source/gInfer/Grammar.cpp
        ${CMAKE_SOURCE_DIR}/source/gInfer/Grammar.h
        ${CMAKE_SOURCE_DIR}/source/gInfer/Rule.cpp
        ${CMAKE_SOURCE_DIR}/source/gInfer/Rule.h
        ${CMAKE_SOURCE_DIR}/source/gInfer/Symbol.cpp
        ${CMAKE_SOURCE_DIR}/source/gInfer/Symbol.h
        pybinding.cpp)

target_link_libraries(pyGInfer
        PRIVATE gInfer
        )