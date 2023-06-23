cmake_minimum_required(VERSION 3.15)

# Dependencies
# ============

if(FS_DOWNLOAD_GBENCHMARK OR FS_GBENCHMARK_SRC_DIR)
    include(FetchContent)

    if(FS_DOWNLOAD_GBENCHMARK)
        message(STATUS "Downloading google-benchmark v1.8.0")

        FetchContent_Declare(googlebenchmark
            GIT_REPOSITORY https://github.com/google/benchmark.git
            GIT_TAG v1.8.0)
    else()
        message(STATUS "Build google-benchmark from local directory")

        FetchContent_Declare(googlebenchmark
            SOURCE_DIR FS_GBENCHMARK_SRC_DIR
            )
    endif()

    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable building benchmark tests" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Disable installing benchmark" FORCE)

    FetchContent_MakeAvailable(googlebenchmark)

else()
    find_package(benchmark REQUIRED)
endif()

# Compiler configuration
# ======================

message(STATUS "Forcing benchmark build type to Release")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
    CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Intel")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-parameter -Wextra -Wreorder")

    # Enable link time optimization and set the default symbol
    # visibility to hidden (very important to obtain small binaries)
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
        # Default symbol visibility
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")

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
  benchmark_basin_graph.cpp
  benchmark_diffusion_adi.cpp
  benchmark_flow_routing.cpp
  benchmark_multi_flow_router.cpp
  benchmark_pflood.cpp
  benchmark_profile_grid.cpp
  benchmark_raster_grid.cpp
  benchmark_single_flow_router.cpp
  benchmark_sink_resolver.cpp
  benchmark_spl.cpp
)

# -- build a target for each benchmark
foreach(filename IN LISTS FASTSCAPELIB_BENCHMARK_SRC)
    string(REPLACE ".cpp" "" targetname ${filename})
    add_executable(${targetname} ${filename} ${FASTSCAPELIB_BENCHMARK_HEADERS} main.cpp)
    target_link_libraries(${targetname} PRIVATE fastscapelib benchmark::benchmark)
    target_compile_features(${targetname} PRIVATE cxx_std_17)
endforeach()

# -- build a global target for all benchmarks
add_executable(benchmark_fastscapelib
        ${FASTSCAPELIB_BENCHMARK_SRC} ${FASTSCAPELIB_BENCHMARK_HEADERS} main.cpp)

target_link_libraries(benchmark_fastscapelib
        PRIVATE fastscapelib benchmark::benchmark)

target_compile_features(benchmark_fastscapelib PRIVATE cxx_std_17)
