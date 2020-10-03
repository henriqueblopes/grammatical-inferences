add_library(source INTERFACE)

target_include_directories(source
        INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/source>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )


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
