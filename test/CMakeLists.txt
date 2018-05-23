cmake_minimum_required(VERSION 3.1)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(fastscapelib-test)

    find_package(xtensor REQUIRED CONFIG)
    set(XTENSOR_INCLUDE_DIR ${xtensor_INCLUDE_DIR})
    find_package(fastscapelib REQUIRED CONFIG)
    set(FASTSCAPELIB_INCLUDE_DIR ${fastscapelib_INCLUDE_DIRS})
endif ()

message(STATUS "Forcing tests build type to Release")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

include(CheckCXXCompilerFlag)

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter -Wextra -Wreorder -Wconversion")

    CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)

    if (HAS_CPP14_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
      # cleaner?
      # see https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right/
      # cmake min version: 3.8
      #target_compile_features(fastscapelib INTERFACE cxx_std_14)
    else()
      message(FATAL_ERROR "Unsupported compiler -- fastscapelib requires C++14 support!")
    endif()
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj")
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif()

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    if(DOWNLOAD_GTEST)
        # Download and unpack googletest at configure time
        configure_file(downloadGTest.cmake.in googletest-download/CMakeLists.txt)
    else()
        # Copy local source of googletest at configure time
        configure_file(copyGTest.cmake.in googletest-download/CMakeLists.txt)
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                     ${CMAKE_CURRENT_BINARY_DIR}/googletest-build)

    set(GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include")
    set(GTEST_BOTH_LIBRARIES  gtest_main gtest)
else()
    find_package(GTest REQUIRED)
endif()

find_package(Threads)

message(STATUS "Prefix path: " ${CMAKE_PREFIX_PATH})

find_package(xtensor REQUIRED)

include_directories(${XTENSOR_INCLUDE_DIR})
include_directories(${FASTSCAPELIB_INCLUDE_DIR})
include_directories(${GTEST_INCLUDE_DIRS})

set(FASTSCAPELIB_TESTS
  main.cpp
  test_flow_routing.cpp
  test_sinks.cpp
)

set(FASTSCAPELIB_TARGET test_fastscapelib)

add_executable(${FASTSCAPELIB_TARGET} ${FASTSCAPELIB_TESTS} ${FASTSCAPELIB_HEADERS} ${XTENSOR_HEADERS})
if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    add_dependencies(${FASTSCAPELIB_TARGET} gtest_main)
endif()

target_link_libraries(${FASTSCAPELIB_TARGET} fastscapelib ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_custom_target(run_tests COMMAND test_fastscapelib DEPENDS ${FASTSCAPELIB_TARGET})
