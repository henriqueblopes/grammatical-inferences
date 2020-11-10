CPMAddPackage(
        NAME Catch2
        GITHUB_REPOSITORY catchorg/Catch2
        VERSION 2.12.1
        OPTIONS
        "CATCH_USE_VALGRIND OFF"
        "CATCH_BUILD_TESTING OFF"
        "CATCH_BUILD_EXAMPLES OFF"
        "CATCH_BUILD_EXTRA_TESTS OFF"
        "CATCH_ENABLE_COVERAGE OFF"
        "CATCH_ENABLE_WERROR ON"
        "CATCH_INSTALL_DOCS OFF"
        "CATCH_INSTALL_HELPERS OFF")
include(${Catch2_SOURCE_DIR}/contrib/Catch.cmake)
add_library(CATCH INTERFACE)

target_include_directories(CATCH
        INTERFACE
        test/
        source/)

add_executable(unit_tests
        main.cpp
        testGrammar.cpp)

target_link_libraries(unit_tests
        Catch2::Catch2
        gInfer
        )

catch_discover_tests(unit_tests)



