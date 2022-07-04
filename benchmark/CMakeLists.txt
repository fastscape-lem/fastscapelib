cmake_minimum_required(VERSION 3.1)

message(STATUS "Forcing benchmark build type to Release")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)


# Dependencies
# ============

# -- Build benchmark independently of fastscapelib
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(xtensor-benchmark)

    find_package(xtensor REQUIRED CONFIG)
    find_package(fastscapelib REQUIRED CONFIG)
endif ()

# -- Google benchmark
if(DOWNLOAD_GBENCHMARK OR GBENCHMARK_SRC_DIR)
    if(DOWNLOAD_GBENCHMARK)
        # Download and unpack googlebenchmark at configure time
        configure_file(downloadGBenchmark.cmake.in
        googlebenchmark-download/CMakeLists.txt)
    else()
        # Copy local source of googlebenchmark at configure time
        configure_file(copyGBenchmark.cmake.in
        googlebenchmark-download/CMakeLists.txt)
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-download)

    if(result)
        message(FATAL_ERROR "CMake step for googlebenchmark failed: ${result}")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-download)

    if(result)
        message(FATAL_ERROR "Build step for googlebenchmark failed: ${result}")
    endif()

    # Add googlebenchmark directly to our build.
    add_subdirectory(
        ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-src
        ${CMAKE_CURRENT_BINARY_DIR}/googlebenchmark-build
        EXCLUDE_FROM_ALL)

    set(GBENCHMARK_LIBRARIES benchmark)

else()
    find_package(benchmark REQUIRED)
endif()

# Compiler configuration
# ======================

include(CheckCXXCompilerFlag)

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
    CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Intel")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-parameter -Wextra -Wreorder")

    CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)

    if (HAS_CPP14_FLAG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
        # cleaner?
        # https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right/
        # cmake min version: 3.8
        #target_compile_features(fastscapelib INTERFACE cxx_std_14)
    else()
        message(FATAL_ERROR "Unsupported compiler -- fastscapelib requires C++14 support!")
    endif()

    # Enable link time optimization and set the default symbol
    # visibility to hidden (very important to obtain small binaries)
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
        # Default symbol visibility
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")

        # Check for Link Time Optimization support
        # (GCC/Clang)
        # LTO had to be removed as google benchmark doesn't build with it
        # CHECK_CXX_COMPILER_FLAG("-flto" HAS_LTO_FLAG)
        # if (HAS_LTO_FLAG)
        #     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
        # endif()

        # Intel equivalent to LTO is called IPO
        if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
            CHECK_CXX_COMPILER_FLAG("-ipo" HAS_IPO_FLAG)
            if(HAS_IPO_FLAG)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ipo")
            endif()
        endif()
    endif()
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj")
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
    foreach(flag_var
            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        string(REPLACE "/MD" "-MT" ${flag_var} "${${flag_var}}")
    endforeach()
endif()

# Build
# =====

set(FASTSCAPELIB_BENCHMARK_HEADERS benchmark_setup.hpp)

set(FASTSCAPELIB_BENCHMARK_SRC
  benchmark_diffusion_adi.cpp
  benchmark_flow_routing.cpp
  benchmark_pflood.cpp
  benchmark_profile_grid.cpp
  benchmark_raster_grid.cpp
  benchmark_single_flow_router.cpp
  benchmark_spl.cpp
)

# -- build a target for each benchmark
foreach(filename IN LISTS FASTSCAPELIB_BENCHMARK_SRC)
    string(REPLACE ".cpp" "" targetname ${filename})
    add_executable(${targetname} ${filename} ${FASTSCAPELIB_BENCHMARK_HEADERS} main.cpp)
    target_link_libraries(${targetname} PRIVATE fastscapelib ${GBENCHMARK_LIBRARIES})
endforeach()

# -- build a global target for all benchmarks
add_executable(benchmark_fastscapelib
        ${FASTSCAPELIB_BENCHMARK_SRC} ${FASTSCAPELIB_BENCHMARK_HEADERS} main.cpp)

target_link_libraries(benchmark_fastscapelib
        PRIVATE fastscapelib ${GBENCHMARK_LIBRARIES})
