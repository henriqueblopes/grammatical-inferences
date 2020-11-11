#[[add_library(source INTERFACE)]]

add_library(gInfer STATIC
        gInfer/Grammar.h
        gInfer/Grammar.cpp
        gInfer/Rule.h
        gInfer/Rule.cpp
        gInfer/Symbol.h
        gInfer/Symbol.cpp
        )

set_target_properties(gInfer PROPERTIES VERSION ${PROJECT_VERSION})

#[[set_target_properties(gInfer PROPERTIES PUBLIC_HEADER Grammar.h Symbol.h Rule.h)]]
set(MY_HEADERS gInfer/Grammar.h gInfer/Rule.h gInfer/Symbol.h)

target_include_directories(gInfer PUBLIC .)

include(GNUInstallDirs)

install(FILES ${MY_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT devel)

install(TARGETS gInfer
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


