cmake_minimum_required(VERSION 3.1)
project(fastscapelib)

set(CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules
  ${CMAKE_MODULE_PATH})

set(FASTSCAPELIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Versioning
# ==========

include(VersionUtils)

get_version_pieces(VERSION_PIECES)

format_version_pep440(VERSION_STR VERSION_PIECES)

get_version_numbers(VERSION_MAJOR VERSION_MINOR VERSION_PATCH
  VERSION_PIECES)

list(GET VERSION_PIECES 0 GIT_HASH_FULL)

message(STATUS "Building ${PROJECT_NAME} version ${VERSION_STR}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/version.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/fastscapelib/version.hpp
  @ONLY)

set(VERSION_HEADER_FILE
  "${CMAKE_CURRENT_BINARY_DIR}/include/fastscapelib/version.hpp")

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

# Dependencies
# ============

# TODO: xtl should be found by xtensor config.
# fix upstream"find_dependency(xtl REQUIRED)" in xtensorConfig.cmake?
find_package(xtl REQUIRED)
message(STATUS "Found xtl: ${xtl_INCLUDE_DIRS}/xtl")

find_package(xtensor REQUIRED)
message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}/xtensor")

# Build
# =====

set(FASTSCAPELIB_HEADERS
  ${VERSION_HEADER_FILE}
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/fastscapelib.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/consts.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/sinks.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow_routing.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/bedrock_channel.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/hillslope.hpp
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
  VERSION ${VERSION_STR}
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
