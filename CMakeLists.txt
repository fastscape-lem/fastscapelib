cmake_minimum_required(VERSION 3.1)
project(fastscapelib)

set(FASTSCAPELIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Versionning
# ===========

file(
  STRINGS
    "${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/config.hpp"
    fastscapelib_version_defines
  REGEX
    "#define FASTSCAPELIB_VERSION_(MAJOR|MINOR|PATCH)"
)

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

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# TODO: xtl should be found by xtensor config.
# fix upstream"find_dependency(xtl REQUIRED)" in xtensorConfig.cmake?
find_package(xtl REQUIRED)
message(STATUS "Found xtl: ${xtl_INCLUDE_DIRS}/xtl")

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
)

set(FASTSCAPELIB_TARGET fastscapelib)

add_library(${FASTSCAPELIB_TARGET} INTERFACE)

target_include_directories(${FASTSCAPELIB_TARGET}
  INTERFACE
    $<BUILD_INTERFACE:${FASTSCAPELIB_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(${FASTSCAPELIB_TARGET} INTERFACE xtensor)

OPTION(BUILD_TESTS "build fastscapelib test suite" OFF)
OPTION(BUILD_BENCHMARK "build fastscapelib benchmark suite" OFF)
OPTION(DOWNLOAD_GTEST "build gtest from downloaded sources" OFF)
OPTION(DOWNLOAD_GBENCHMARK "download google benchmark and build from source" ON)
OPTION(BUILD_PYTHON_MODULE "build fastscapelib Python bindings" OFF)

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
  set(BUILD_TESTS ON)
endif()

if(BUILD_TESTS)
  include_directories(${FASTSCAPELIB_INCLUDE_DIR})
  add_subdirectory(test)
endif()

if(BUILD_BENCHMARK)
    add_subdirectory(benchmark)
endif()

if(BUILD_PYTHON_MODULE)
  include_directories(${FASTSCAPELIB_INCLUDE_DIR})
  add_subdirectory(python)
endif()


# Installation
# ============

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# -- install header files
install(FILES ${FASTSCAPELIB_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fastscapelib
)

# -- create/install targets (also importable from the build dir)
install(TARGETS ${FASTSCAPELIB_TARGET}
  EXPORT ${PROJECT_NAME}-targets
)

install(EXPORT ${PROJECT_NAME}-targets
  FILE ${PROJECT_NAME}Targets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

export(EXPORT ${PROJECT_NAME}-targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake
)

# -- create/install cmake config
configure_package_config_file(${PROJECT_NAME}Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# fastscapelib is header-only and does not depend on the architecture.
# Remove CMAKE_SIZEOF_VOID_P from fastscapelibConfigVersion.cmake
# so that a fastscapelibConfig.cmake generated for a 64 bit target
# can be used for 32 bit targets and vice versa.
set(_FASTSCAPELIB_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY AnyNewerVersion
)
set(CMAKE_SIZEOF_VOID_P ${_FASTSCAPELIB_CMAKE_SIZEOF_VOID_P})

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# -- create/install pkg-config
configure_file(${PROJECT_NAME}.pc.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  @ONLY
)
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/"
)
