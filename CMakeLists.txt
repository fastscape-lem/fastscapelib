cmake_minimum_required(VERSION 3.15)

project(
  fastscapelib
  VERSION "0.2.1"
  LANGUAGES CXX)

message(STATUS "Building fastscapelib v${PROJECT_VERSION}")

set(FASTSCAPELIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

OPTION(FS_DOWNLOAD_XTENSOR "Download and use xtensor development version" OFF)
OPTION(FS_DOWNLOAD_XTENSOR_PYTHON "Download xtensor-python (python bindings)" OFF)
OPTION(FS_BUILD_TESTS "Build fastscapelib test suite" OFF)
OPTION(FS_DOWNLOAD_GTEST "Build gtest from downloaded sources" OFF)
OPTION(FS_GTEST_SRC_DIR "Build gtest from local sources" OFF)
OPTION(FS_BUILD_BENCHMARKS "Build fastscapelib benchmark suite" OFF)
OPTION(FS_DOWNLOAD_GBENCHMARK "Build gbenchmark from downloaded sources" OFF)
OPTION(FS_GBENCHMARK_SRC_DIR "Build gbenchmark from local sources" OFF)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH})

# Dependencies
# ============

include(FetchContent)

# -- get the whole xtensor stack (only when building python bindings)
if(SKBUILD AND FS_DOWNLOAD_XTENSOR_PYTHON)
    message(STATUS "Downloading xtl, xtensor and xtensor-python libraries")

    FetchContent_Declare(xtl
        GIT_REPOSITORY https://github.com/xtensor-stack/xtl
        GIT_TAG 0.7.5
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(xtl)

    FetchContent_Declare(xtensor
        GIT_REPOSITORY https://github.com/xtensor-stack/xtensor
        GIT_TAG 0.24.6
        GIT_SHALLOW TRUE)
    set(CPP17 ON CACHE BOOL "Enable C++17 for xtensor" FORCE)
    FetchContent_MakeAvailable(xtensor)

    FetchContent_Declare(xtensor-python
        GIT_REPOSITORY https://github.com/xtensor-stack/xtensor-python
        GIT_TAG 0.26.1
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(xtensor-python)

    set_property(DIRECTORY ${xtl_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
    set_property(DIRECTORY ${xtensor_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
    set_property(DIRECTORY ${xtensor-python_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)

# -- get xtensor development version (need xtl already installed)
elseif(FS_DOWNLOAD_XTENSOR)
    message(STATUS "Downloading xtensor development version (master branch)")

    FetchContent_Declare(xtensor
      GIT_REPOSITORY https://github.com/xtensor-stack/xtensor
      GIT_TAG master
      GIT_SHALLOW TRUE)

    set(CPP17 ON CACHE BOOL "Enable C++17 for xtensor" FORCE)
    set(DOWNLOAD_GBENCHMARK OFF CACHE BOOL
      "skip downloading google-benchmark via xtensor" FORCE)

    FetchContent_MakeAvailable(xtensor)

else()
    find_package(xtensor REQUIRED)
    message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}/xtensor")
endif()

# Installation directories
# ========================

include(GNUInstallDirs)
message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})

set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for fastscapelibConfig.cmake")

message("  -> binaries: " ${CMAKE_INSTALL_BINDIR})
message("  -> libs: " ${CMAKE_INSTALL_LIBDIR})
message("  -> includes: " ${CMAKE_INSTALL_INCLUDEDIR})
message("  -> cmake config: " ${CMAKECONFIG_INSTALL_DIR})

include(CMakePackageConfigHelpers)

# Configure/Build
# ===============

# -- fastscapelib
set(FASTSCAPELIB_HEADERS
  # Algo
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/algo/pflood.hpp
  # Eroders
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/eroders/diffusion_adi.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/eroders/spl.hpp
  # Flow
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow/basin_graph.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow/flow_graph_impl.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow/flow_graph.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow/flow_operator.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow/flow_router.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow/flow_snapshot.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow/sink_resolver.hpp
  # Grid
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/grid/base.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/grid/profile_grid.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/grid/raster_grid.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/grid/structured_grid.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/grid/trimesh.hpp
  # Utils
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils/consts.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils/iterators.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils/union_find.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils/utils.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils/xtensor_utils.hpp
  # Version
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/version.hpp
)

set(FASTSCAPELIB_TARGET fastscapelib)

add_library(${FASTSCAPELIB_TARGET} INTERFACE)

target_include_directories(${FASTSCAPELIB_TARGET}
  INTERFACE
    $<BUILD_INTERFACE:${FASTSCAPELIB_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(${FASTSCAPELIB_TARGET} INTERFACE xtensor)

target_compile_features(${FASTSCAPELIB_TARGET} INTERFACE cxx_std_17)

# -- optional subdirectories
if(FS_DOWNLOAD_GTEST OR FS_GTEST_SRC_DIR)
  set(FS_BUILD_TESTS ON)
endif()

if(FS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(FS_DOWNLOAD_GBENCHMARK OR FS_GBENCHMARK_SRC_DIR)
  set(FS_BUILD_BENCHMARKS ON)
endif()

if(FS_BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif()

if(SKBUILD)
  add_subdirectory(python)
endif()

# -- CMake configuration
#    o Package config (makes the package importable using find_package in another project)
configure_package_config_file(cmake/${PROJECT_NAME}Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})

# Installation of the C++ library
# ===============================

if(NOT SKBUILD)
  # -- fastscapelib header files
  install(FILES ${FASTSCAPELIB_HEADERS}
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fastscapelib)

  # -- CMake configuration
  #    o create/install targets
  install(TARGETS ${FASTSCAPELIB_TARGET}
          EXPORT ${PROJECT_NAME}-targets)

  install(EXPORT ${PROJECT_NAME}-targets
          FILE ${PROJECT_NAME}Targets.cmake
          DESTINATION ${CMAKECONFIG_INSTALL_DIR})

  #    o make targets also importable from the build dir
  export(EXPORT ${PROJECT_NAME}-targets
         FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)

  #    o version file (fastscapelib is header-only and does not depend on the architecture)
  set(_FASTSCAPELIB_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})
  unset(CMAKE_SIZEOF_VOID_P)
  write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
                                   VERSION ${PROJECT_VERSION}
                                   COMPATIBILITY AnyNewerVersion)
  set(CMAKE_SIZEOF_VOID_P ${_FASTSCAPELIB_CMAKE_SIZEOF_VOID_P})

  #    o install config files
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()
