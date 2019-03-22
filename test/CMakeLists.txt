cmake_minimum_required(VERSION 3.1)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  project(fastscapelib-test)

  find_package(xtensor REQUIRED CONFIG)
  set(XTENSOR_INCLUDE_DIR ${xtensor_INCLUDE_DIR})
  find_package(fastscapelib REQUIRED CONFIG)
  set(FASTSCAPELIB_INCLUDE_DIR ${fastscapelib_INCLUDE_DIRS})
endif ()

# message(STATUS "Forcing tests build type to Release")
# set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

include(CheckCXXCompilerFlag)

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
    CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
    (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND NOT WIN32))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wreorder")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Wsign-conversion")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast -Wunused-variable")

  CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)

  if (HAS_CPP14_FLAG)
    set(CMAKE_CXX_STANDARD 14)
    message(STATUS "Building with -std=c++14")
    # cleaner?
    # https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right/
    # cmake min version: 3.8
    #target_compile_features(fastscapelib INTERFACE cxx_std_14)
  else()
    message(FATAL_ERROR "Unsupported compiler -- fastscapelib requires C++14 support!")
  endif()

  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
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
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
    )

  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
    )

  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()

  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(
    ${CMAKE_CURRENT_BINARY_DIR}/googletest-src
    ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
    EXCLUDE_FROM_ALL
    )

  set(GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include")
  set(GTEST_BOTH_LIBRARIES gtest_main gtest)

else()
  find_package(GTest REQUIRED)
endif()

include_directories(${GTEST_INCLUDE_DIRS} SYSTEM)

set(FASTSCAPELIB_TEST_SRC
  main.cpp
  test_grid.cpp
  test_hillslope.cpp
  test_flow_routing.cpp
  test_sinks.cpp
  test_bedrock_channel.cpp
  )

set(CMAKE_VERBOSE_MAKEFILE ON)

set(FASTSCAPELIB_TEST_TARGET test_fastscapelib)

add_executable(${FASTSCAPELIB_TEST_TARGET} ${FASTSCAPELIB_TEST_SRC})

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    add_dependencies(${FASTSCAPELIB_TEST_TARGET} gtest_main)
endif()

target_include_directories(${FASTSCAPELIB_TEST_TARGET}
  PRIVATE ${FASTSCAPELIB_INCLUDE_DIR})

target_link_libraries(${FASTSCAPELIB_TEST_TARGET}
  fastscapelib ${GTEST_BOTH_LIBRARIES})

add_custom_target(run_tests
  ALL
  COMMAND ${FASTSCAPELIB_TEST_TARGET}
  DEPENDS ${FASTSCAPELIB_TEST_TARGET}
)
