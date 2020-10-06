cmake_minimum_required(VERSION 3.14)
project(grammatical-inferences)

set(CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")

option(BUILD_WITH_PEDANTIC_WARNINGS "Use pedantic warnings. This is useful for developers because many of these warnings will be in continuous integration anyway." ${DEBUG_MODE})
option(BUILD_INSTALLER "Build an installation package" ${MASTER_PROJECT})
option(BUILD_PACKAGE "Build an installation package" ${MASTER_PROJECT})

include_directories(${PYTHON_INCLUDE_DIRS})

#[[add_subdirectory(source bindings)]]

add_executable(grammatical-inferences
        source/Grammar.cpp
        source/Grammar.h
        source/Rule.cpp
        source/Rule.h
        source/Symbol.cpp
        source/Symbol.h
        source/main.cpp
        source/InputWords.cpp
        source/InputWords.h)


#[[
find_package(pybind11 REQUIRED)

pybind11_add_module(klpcsg
        source/Grammar.cpp
        source/Grammar.h
        source/Rule.cpp
        source/Rule.h
        source/Symbol.cpp
        source/Symbol.h
        pybinding.cpp)
]]


target_link_libraries(grammatical-inferences stdc++fs ${PYTHON_LIBRARIES})# pybind11::embed)

if (BUILD_INSTALLER)
    set(CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/source)

    write_basic_package_version_file(
            sourceConfigVersion.cmake
            VERSION ${PACKAGE_VERSION}
            COMPATIBILITY AnyNewerVersion
    )
    configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/sourceConfig.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/sourceConfig.cmake
            INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/source)

    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/sourceConfig.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/sourceConfigVersion.cmake
            COMPONENT "CPP_Library"
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/source
            )
endif()


if (BUILD_INSTALLER AND BUILD_PACKAGE)
    set(CPACK_PACKAGE_NAME "grammatical-inferences")
    set(CPACK_PACKAGE_VENDOR "https://gitlab.com/alandefreitas")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Grammatical Inferences for for C++/Python")
    set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
    set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
    set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

    FILE(READ .gitignore GITIGNORE_CONTENTS)
    STRING(REGEX REPLACE ";" "\\\\;" GITIGNORE_CONTENTS "${GITIGNORE_CONTENTS}")
    STRING(REGEX REPLACE "\n" ";" GITIGNORE_CONTENTS "${GITIGNORE_CONTENTS}")
    set(CPACK_SOURCE_IGNORE_FILES ${GITIGNORE_CONTENTS})

    include(CPack)
endif()