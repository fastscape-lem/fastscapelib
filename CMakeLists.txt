cmake_minimum_required(VERSION 3.1)
project(fastscapelib)

set(CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules
  ${CMAKE_MODULE_PATH})

set(FASTSCAPELIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Versioning
# ==========

include(GetGitRevisionDescription)

get_git_head_revision(_ GIT_HASH_FULL)
git_describe(GIT_REV_DESCRIPTION --tags --always --dirty)

if(${GIT_REV_DESCRIPTION} MATCHES
    "^[v]*([0-9\\.]+)\\-([0-9]+)\\-g([A-Za-z0-9]+)([dirty\\-]*).*")
  # found tag + additional commits
  set(GIT_CLOSEST_TAG ${CMAKE_MATCH_1})
  set(GIT_COMMIT_COUNT ${CMAKE_MATCH_2})
  set(GIT_HASH_SHORT ${CMAKE_MATCH_3})
  set(GIT_LOCAL_CHANGES ${CMAKE_MATCH_4})

elseif(${GIT_REV_DESCRIPTION} MATCHES
    "^[v]*([0-9]+\\.[0-9]+\\.[0-9]+)([dirty\\-]*).*")
  # HEAD is a tag
  set(GIT_CLOSEST_TAG ${CMAKE_MATCH_1})
  set(GIT_COMMIT_COUNT "")
  set(GIT_HASH_SHORT "")
  set(GIT_LOCAL_CHANGES ${CMAKE_MATCH_2})

elseif(${GIT_REV_DESCRIPTION} MATCHES
    "^([A-Za-z0-9]+)([dirty\\-]*).*")
  # no tag found
  set(GIT_CLOSEST_TAG "0+untagged")
  git_get_commit_count(GIT_COMMIT_COUNT)
  set(GIT_HASH_SHORT ${CMAKE_MATCH_1})
  set(GIT_LOCAL_CHANGES ${CMAKE_MATCH_2})

else()
  message(FATAL_ERROR
    "Could not parse version correctly. Output of git describe: "
    ${GIT_REV_DESCRIPTION})

endif()

if(${GIT_CLOSEST_TAG} MATCHES "([0-9]+)\\.([0-9]+)\\.([0-9]+)")
  set(VERSION_MAJOR ${CMAKE_MATCH_1})
  set(VERSION_MINOR ${CMAKE_MATCH_2})
  set(VERSION_PATCH ${CMAKE_MATCH_3})
else()
  set(VERSION_MAJOR "0")
  set(VERSION_MINOR "0")
  set(VERSION_PATCH "0")
endif()

if(NOT GIT_HASH_SHORT)
  set(VERSION_GIT_INFO "")
else()
  set(VERSION_GIT_INFO +${GIT_COMMIT_COUNT}.g${GIT_HASH_SHORT})
endif()

set(VERSION_STR ${GIT_CLOSEST_TAG}${VERSION_GIT_INFO}${GIT_LOCAL_CHANGES})
message(STATUS "Building ${PROJECT_NAME} version ${VERSION_STR}")

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
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/fastscapelib.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/config.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/consts.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/sinks.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow_routing.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/bedrock_channel.hpp
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
