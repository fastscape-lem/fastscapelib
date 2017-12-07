cmake_minimum_required(VERSION 3.1)
project(fastscape)

set(FASTSCAPE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Versionning
# ===========

file(STRINGS "${FASTSCAPE_INCLUDE_DIR}/fastscape/config.hpp" fastscape_version_defines
     REGEX "#define FASTSCAPE_VERSION_(MAJOR|MINOR|PATCH)")
foreach(ver ${fastscape_version_defines})
    if(ver MATCHES "#define FASTSCAPE_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(FASTSCAPE_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()
set(${PROJECT_NAME}_VERSION
    ${FASTSCAPE_VERSION_MAJOR}.${FASTSCAPE_VERSION_MINOR}.${FASTSCAPE_VERSION_PATCH})
message(STATUS "Building fastscape v${${PROJECT_NAME}_VERSION}")

# Dependencies
# ============

find_package(xtensor REQUIRED)
message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}/xtensor")

# Build
# =====

set(FASTSCAPE_HEADERS
  ${FASTSCAPE_INCLUDE_DIR}/fastscape/config.hpp
  ${FASTSCAPE_INCLUDE_DIR}/fastscape/utils.hpp
  ${FASTSCAPE_INCLUDE_DIR}/fastscape/consts.hpp
  ${FASTSCAPE_INCLUDE_DIR}/fastscape/sinks.hpp
)

include_directories(${xtensor_INCLUDE_DIRS})

OPTION(BUILD_TESTS "fastscape test suite" OFF)
OPTION(DOWNLOAD_GTEST "build gtest from downloaded sources" OFF)

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    set(BUILD_TESTS ON)
endif()

if(BUILD_TESTS)
    include_directories(${FASTSCAPE_INCLUDE_DIR})
    add_subdirectory(test)
endif()

# Installation
# ============

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(FILES ${FASTSCAPE_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fastscape)

set(FASTSCAPE_CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for fastscapeConfig.cmake")

configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION ${FASTSCAPE_CMAKECONFIG_INSTALL_DIR})

# fastscape is header-only and does not depend on the architecture.
# Remove CMAKE_SIZEOF_VOID_P from fastscapeConfigVersion.cmake so that an fastscapeConfig.cmake
# generated for a 64 bit target can be used for 32 bit targets and vice versa.
set(_FASTSCAPE_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
                                 VERSION ${${PROJECT_NAME}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
set(CMAKE_SIZEOF_VOID_P ${_FASTSCAPE_CMAKE_SIZEOF_VOID_P})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${FASTSCAPE_CMAKECONFIG_INSTALL_DIR})

configure_file(${PROJECT_NAME}.pc.in
               "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
                @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")
