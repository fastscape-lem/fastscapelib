cmake_minimum_required(VERSION 3.1)

function(set_version_str INCLUDEDIR)
    file(STRINGS "${INCLUDEDIR}/fastscapelib/version.hpp" fscape_version_defines
        REGEX "#define FASTSCAPELIB_VERSION_(MAJOR|MINOR|PATCH)")
    foreach(ver ${fscape_version_defines})
        if(ver MATCHES "#define FASTSCAPELIB_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
            set(FASTSCAPELIB_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
        endif()
    endforeach()
    set(VERSION_STR
        ${FASTSCAPELIB_VERSION_MAJOR}.${FASTSCAPELIB_VERSION_MINOR}.${FASTSCAPELIB_VERSION_PATCH} PARENT_SCOPE)
endfunction()
