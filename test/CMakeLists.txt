#[[CPMAddPackage(
        NAME Catch2
        GITHUB_REPOSITORY catchorg/Catch2
        VERSION 3.2.1
        OPTIONS
        "CATCH_USE_VALGRIND OFF"
        "CATCH_BUILD_TESTING OFF"
        "CATCH_BUILD_EXAMPLES OFF"
        "CATCH_BUILD_EXTRA_TESTS OFF"
        "CATCH_ENABLE_COVERAGE OFF"
        "CATCH_ENABLE_WERROR ON"
        "CATCH_INSTALL_DOCS OFF"
        "CATCH_INSTALL_HELPERS OFF")
include(${Catch2_SOURCE_DIR}/contrib/Catch.cmake)]]

Include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.0.1 # or a later release
)

FetchContent_MakeAvailable(Catch2)

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
include(Catch)
catch_discover_tests(unit_tests)



