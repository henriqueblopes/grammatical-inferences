Include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.13.1)
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

include(CTest)
