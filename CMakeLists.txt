cmake_minimum_required(VERSION 3.1)

project(fastscapelib)

set(FASTSCAPELIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

OPTION(DOWNLOAD_XTENSOR "Download xtensor from github (master)" OFF)
OPTION(DOWNLOAD_GTEST "build gtest from downloaded sources" OFF)
OPTION(BUILD_TESTS "build fastscapelib test suite" OFF)
OPTION(BUILD_BENCHMARK "build fastscapelib benchmark suite" OFF)
OPTION(BUILD_PYTHON_MODULE "build fastscapelib Python bindings" OFF)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH})

# Versioning
# ==========

include(VersionUtils)
set_version_str("${FASTSCAPELIB_INCLUDE_DIR}")
message(STATUS "Building fastscapelib v${VERSION_STR}")

# Dependencies
# ============

# -- xtensor (eventually from source download)
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
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/version.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/fastscapelib.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/utils.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/profile_grid.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/raster_grid.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/consts.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/sinks.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow_routing.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow_graph.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow_router.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/flow_router_factory.hpp
  ${FASTSCAPELIB_INCLUDE_DIR}/fastscapelib/sink_resolver.hpp
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

target_compile_features(${FASTSCAPELIB_TARGET} INTERFACE cxx_std_14)

# -- optional subdirectories
if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
  set(BUILD_TESTS ON)
endif()

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(BUILD_BENCHMARK)
  add_subdirectory(benchmark)
endif()

if(BUILD_PYTHON_MODULE)
  add_subdirectory(python)
endif()

# -- CMake configuration

#    o Package config (makes the package importable using find_package in another project)
configure_package_config_file(${PROJECT_NAME}Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})

#    o pkg-config
configure_file(${PROJECT_NAME}.pc.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
    @ONLY)

# Installation
# ============

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
                                 VERSION ${VERSION_STR}
                                 COMPATIBILITY AnyNewerVersion)
set(CMAKE_SIZEOF_VOID_P ${_FASTSCAPELIB_CMAKE_SIZEOF_VOID_P})

#    o Package config
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION ${CMAKECONFIG_INSTALL_DIR})

#    o pkg-config
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")
