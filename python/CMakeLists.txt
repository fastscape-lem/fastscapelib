cmake_minimum_required(VERSION 3.1)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  project(fastscapelib-python)

  find_package(xtensor REQUIRED CONFIG)
  find_package(fastscapelib REQUIRED CONFIG)
endif ()

# Dependencies
# ============

find_package(xtensor-python REQUIRED)
message(STATUS "Found xtensor-python: ${xtensor-python_INCLUDE_DIRS}/xtensor-python")

find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

find_package(NumPy REQUIRED)
message(STATUS "Found NumPy: v${NUMPY_VERSION} ${NUMPY_INCLUDE_DIRS}")

# Build
# =====

# TODO: (upstream) modern cmake (target/properties) for those dependencies
include_directories(${xtensor-python_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIRS})

set(FASTSCAPELIB_PYTHON_TARGET _fastscapelib_py)

pybind11_add_module(${FASTSCAPELIB_PYTHON_TARGET}
  src/grid.cpp
  src/main.cpp
  )

target_link_libraries(${FASTSCAPELIB_PYTHON_TARGET}
  INTERFACE xtensor xtensor-python fastscapelib)
