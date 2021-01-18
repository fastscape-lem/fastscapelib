cmake_minimum_required(VERSION 3.1)


# Dependencies
# ============

# -- Build tests independently of fastscapelib
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(fastscapelib-test)

    find_package(xtensor REQUIRED CONFIG)
    find_package(fastscapelib REQUIRED CONFIG)
endif()

# -- GTest (eventually from source dir or download)
if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    if(DOWNLOAD_GTEST) # Download and unpack googletest at configure time
        configure_file(downloadGTest.cmake.in googletest-download/CMakeLists.txt)
    else() # Copy local source of googletest at configure time
        configure_file(copyGTest.cmake.in googletest-download/CMakeLists.txt)
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)

    if(result)
        message(FATAL_ERROR "CMake configuration step for googletest failed: ${result}")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)

    if(result)
        message(FATAL_ERROR "CMake build step for googletest failed: ${result}")
    endif()

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    
    # Add googletest directly to our build. This defines the gtest and gtest_main targets.
    add_subdirectory(
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL
        )

    set(GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include")
    add_library(GTest::GTest INTERFACE IMPORTED)
    target_link_libraries(GTest::GTest INTERFACE gtest)
    add_library(GTest::Main INTERFACE IMPORTED)
    target_link_libraries(GTest::Main INTERFACE gtest_main)
else()
    find_package(GTest REQUIRED)
endif()

# Compiler configuration
# ======================

# -- Set compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
    CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
    (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND NOT WIN32))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wreorder")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Wsign-conversion")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast -Wunused-variable")

    # Check C++14 support (not working with MSVC)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)

    if(NOT HAS_CPP14_FLAG)
        message(FATAL_ERROR "Unsupported compiler -- fastscapelib requires C++14 support!")
    endif()
endif()

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj")
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif()

# Build
# =====

set(FASTSCAPELIB_TEST_SRC
    test_hillslope.cpp
    test_flow_graph.cpp
    test_flow_router.cpp
    test_single_flow_router.cpp
    test_flow_routing.cpp
    test_profile_grid.cpp
    test_raster_grid.cpp
    test_raster_rook_grid.cpp
    test_raster_queen_grid.cpp
    test_raster_bishop_grid.cpp    
    test_sinks.cpp
    test_sink_resolver.cpp
    test_bedrock_channel.cpp
)

set(FASTSCAPELIB_TEST_TARGET test_fastscapelib)

# -- build a target for each benchmark
foreach(filename IN LISTS FASTSCAPELIB_TEST_SRC)
    string(REPLACE ".cpp" "" targetname ${filename})
    add_executable(${targetname} ${filename} main.cpp)
    target_link_libraries(${targetname} PRIVATE fastscapelib GTest::GTest GTest::Main)
endforeach()

# -- build a global target for all benchmarks
add_executable(${FASTSCAPELIB_TEST_TARGET} ${FASTSCAPELIB_TEST_SRC} main.cpp)

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    add_dependencies(${FASTSCAPELIB_TEST_TARGET} gtest_main)
endif()

target_link_libraries(${FASTSCAPELIB_TEST_TARGET} PRIVATE
                      fastscapelib GTest::GTest GTest::Main)

target_compile_features(${FASTSCAPELIB_TEST_TARGET} INTERFACE cxx_std_14)


