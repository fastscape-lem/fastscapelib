cmake_minimum_required(VERSION 3.15)

project(
  fastscapelib
  VERSION "0.1.3"
  LANGUAGES CXX)

message(STATUS "Building fastscapelib v${PROJECT_VERSION}")

set(FASTSCAPELIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

OPTION(FS_DOWNLOAD_XTENSOR "Download and use xtensor development version" OFF)
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

# -- xtensor (optionally from source download)
if(DOWNLOAD_XTENSOR)
    configure_file(downloadXTensor.cmake.in
      xtensor-download/CMakeLists.txt)

    execute_process(
        COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/xtensor-download
    )

    if(result)
        message(FATAL_ERROR "CMake configuration step for xtensor failed: ${result}")
    endif()

    execute_process(
        COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/xtensor-download
    )

    if(result)
        message(FATAL_ERROR "CMake build step for xtensor failed: ${result}")
    endif()

    set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_BINARY_DIR}/xtensor-build ${CMAKE_MODULE_PATH}")
    find_package(xtensor REQUIRED PATHS "${CMAKE_CURRENT_BINARY_DIR}/xtensor-build" NO_DEFAULT_PATH)
else()
    find_package(xtensor REQUIRED)
endif()

get_target_property(xtensor_INCLUDE_DIRS xtensor INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}/xtensor")

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
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/algo/flow_routing.hpp
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
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/grid/unstructured_mesh.hpp
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
configure_package_config_file(${PROJECT_NAME}Config.cmake.in
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
