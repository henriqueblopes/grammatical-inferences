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

#[[target_include_directories(source
        INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/source>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )]]


#[[
if (BUILD_WITH_PEDANTIC_WARNINGS)
    if (MSVC)
        target_compile_options(source INTERFACE /W4 /WX)
    else()
        target_compile_options(source INTERFACE -Wall -Wextra -pedantic -Werror -Wno-self-assign-overloaded)
    endif()
endif()


if(BUILD_INSTALLER)
    install(TARGETS source
            EXPORT sourceTargets
            COMPONENT "CPP_Library"
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/source/
            COMPONENT "CPP_Library"
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            FILES_MATCHING PATTERN "*.h"
            )
    install(EXPORT sourceTargets
            FILE sourceTargets.cmake
            NAMESPACE sourceTargets::
            COMPONENT "CPP_Library"
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/source
            )
endif()
]]
