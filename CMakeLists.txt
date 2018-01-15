cmake_minimum_required(VERSION 3.1)
project(fastscapelib)

set(FASTSCAPELIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Versionning
# ===========

file(STRINGS "${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/config.hpp" fastscapelib_version_defines
     REGEX "#define FASTSCAPELIB_VERSION_(MAJOR|MINOR|PATCH)")
foreach(ver ${fastscapelib_version_defines})
    if(ver MATCHES "#define FASTSCAPELIB_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(FASTSCAPELIB_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()
set(${PROJECT_NAME}_VERSION
    ${FASTSCAPELIB_VERSION_MAJOR}.${FASTSCAPELIB_VERSION_MINOR}.${FASTSCAPELIB_VERSION_PATCH})
message(STATUS "Building fastscapelib v${${PROJECT_NAME}_VERSION}")

# Dependencies
# ============

find_package(xtensor REQUIRED)
message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}/xtensor")

# Build
# =====

set(FASTSCAPELIB_HEADERS
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/fastscapelib.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/config.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/consts.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/sinks.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow_routing.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/basin_graph.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/union_find.hpp
)

include_directories(${xtensor_INCLUDE_DIRS})

OPTION(BUILD_TESTS "fastscapelib test suite" OFF)
OPTION(DOWNLOAD_GTEST "build gtest from downloaded sources" OFF)

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    set(BUILD_TESTS ON)
endif()

if(BUILD_TESTS)
    include_directories(${FASTSCAPELIB_INCLUDE_DIR})
    add_subdirectory(test)
endif()

# Installation
# ============

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(FILES ${FASTSCAPELIB_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fastscapelib)

set(FASTSCAPELIB_CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for fastscapelibConfig.cmake")

configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION ${FASTSCAPELIB_CMAKECONFIG_INSTALL_DIR})

# fastscapelib is header-only and does not depend on the architecture.
# Remove CMAKE_SIZEOF_VOID_P from fastscapelibConfigVersion.cmake so that an fastscapelibConfig.cmake
# generated for a 64 bit target can be used for 32 bit targets and vice versa.
set(_FASTSCAPELIB_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
                                 VERSION ${${PROJECT_NAME}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
set(CMAKE_SIZEOF_VOID_P ${_FASTSCAPELIB_CMAKE_SIZEOF_VOID_P})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${FASTSCAPELIB_CMAKECONFIG_INSTALL_DIR})

configure_file(${PROJECT_NAME}.pc.in
               "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
                @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")
