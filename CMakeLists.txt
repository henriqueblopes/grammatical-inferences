cmake_minimum_required(VERSION 3.14)
project(grammatical-inferences VERSION 1.1.0)

set(CMAKE_CXX_STANDARD 17)


set(CPM_DOWNLOAD_VERSION 0.27.2)
set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")

if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
    message(STATUS "Downloading CPM.cmake")
    file(DOWNLOAD https://github.com/TheLartians/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake ${CPM_DOWNLOAD_LOCATION})
endif()

include(${CPM_DOWNLOAD_LOCATION})
option(CPM_USE_LOCAL_PACKAGES "Try `find_package` before downloading dependencies" ON)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc")

option(BUILD_WITH_PEDANTIC_WARNINGS "Use pedantic warnings. This is useful for developers because many of these warnings will be in continuous integration anyway." OFF)
option(BUILD_INSTALLER "Build an installation package" ON)
option(BUILD_PACKAGE "Build an installation package" ON)
option(BUILD_TESTS "Compile the tests" ON)
option(BUILD_PYTHON_BINDING "Compile the Python binding" ON)

if (BUILD_WITH_PEDANTIC_WARNINGS)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror -Wno-self-assign-overloaded")
endif()



CPMAddPackage(
        NAME pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11/
        VERSION 2.5.0
        OPTIONS "PYBIND11_TEST OFF" "PYBIND11_INSTALL OFF"
)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")

include_directories(${PYTHON_INCLUDE_DIRS})

add_subdirectory(source)

if (BUILD_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(test)
endif()

if (BUILD_PYTHON_BINDING)
    add_subdirectory(bindings)
endif()

find_package(Filesystem REQUIRED)

add_executable(grammatical-inferences
        source/gInfer/Grammar.cpp
        source/gInfer/Grammar.h
        source/gInfer/Rule.cpp
        source/gInfer/Rule.h
        source/gInfer/Symbol.cpp
        source/gInfer/Symbol.h
        examples/main.cpp
        examples/InputWords.cpp
        examples/InputWords.h)


target_link_libraries(grammatical-inferences ${PYTHON_LIBRARIES} pybind11::embed Catch2::Catch2 gInfer std::filesystem)
#[[target_compile_features(grammatical-inferences PRIVATE cxx_std_17)]]


#[[include(ParseAndAddCatchTests)
ParseAndAddCatchTests(unit_tests)]]





if (BUILD_INSTALLER)
    set(CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/GrammaticalInferences)
    message("CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")
    message("CMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}")

    # Create Matplot++ConfigVersion.cmake and install it
    write_basic_package_version_file(
            GrammaticalInferencesConfigVersion.cmake
            VERSION ${PACKAGE_VERSION}
            COMPATIBILITY AnyNewerVersion
    )

    configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/GrammaticalInferences-config.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/GrammaticalInferences-config.cmake
            INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GrammaticalInferences)

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GrammaticalInferences-config.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GrammaticalInferences)

    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindFilesystem.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GrammaticalInferences)

    #[[
    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/sourceConfig.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/sourceConfigVersion.cmake
            COMPONENT "CPP_Library"
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/source/gInfer
            )]]
endif()


if (BUILD_INSTALLER AND BUILD_PACKAGE)
    set(CPACK_PACKAGE_NAME "grammatical-inferences")
    set(CPACK_PACKAGE_VENDOR "https://github.com/henriqueblopes")
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