cmake_minimum_required(VERSION 3.14)
project(grammatical-inferences VERSION 1.1.0)

set(CMAKE_CXX_STANDARD 17)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror -Wno-self-assign-overloaded")


set(CPM_DOWNLOAD_VERSION 0.27.2)
set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")

if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
    message(STATUS "Downloading CPM.cmake")
    file(DOWNLOAD https://github.com/TheLartians/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake ${CPM_DOWNLOAD_LOCATION})
endif()

include(${CPM_DOWNLOAD_LOCATION})
option(CPM_USE_LOCAL_PACKAGES "Try `find_package` before downloading dependencies" ON)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")

option(BUILD_WITH_PEDANTIC_WARNINGS "Use pedantic warnings. This is useful for developers because many of these warnings will be in continuous integration anyway." ${DEBUG_MODE})
option(BUILD_INSTALLER "Build an installation package" ${MASTER_PROJECT})
option(BUILD_PACKAGE "Build an installation package" ${MASTER_PROJECT})



CPMAddPackage(
        NAME pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11/
        VERSION 2.5.0
        OPTIONS "PYBIND11_TEST OFF" "PYBIND11_INSTALL OFF"
)

include_directories(${PYTHON_INCLUDE_DIRS})

Include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.13.1)
FetchContent_MakeAvailable(Catch2)

#[[include(Catch)]]

add_subdirectory(source)

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

find_package(pybind11 REQUIRED)


pybind11_add_module(pyGInfer
        source/Grammar.cpp
        source/Grammar.h
        source/Rule.cpp
        source/Rule.h
        source/Symbol.cpp
        source/Symbol.h
        pybinding.cpp)



target_link_libraries(grammatical-inferences stdc++fs  ${PYTHON_LIBRARIES} pybind11::embed Catch2::Catch2)


add_library(CATCH INTERFACE)
target_include_directories(CATCH
        INTERFACE
            test/)

add_executable(unit_tests
        test/main.cpp
        test/testGrammar.cpp)
target_link_libraries(unit_tests
        Catch2::Catch2
        gInfer
        )

include(CTest)
#[[include(ParseAndAddCatchTests)
ParseAndAddCatchTests(unit_tests)]]

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