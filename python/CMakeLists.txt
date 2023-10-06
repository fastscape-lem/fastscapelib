cmake_minimum_required(VERSION 3.15)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH})

# Dependencies
# ============

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)
message(STATUS "Found python v${Python_VERSION}: ${Python_EXECUTABLE}")
message(STATUS "Found numpy v${Python_NumPy_VERSION}: ${Python_NumPy_INCLUDE_DIRS}")

if(NOT FS_DOWNLOAD_XTENSOR_PYTHON)
  find_package(xtensor-python REQUIRED)
  message(STATUS "Found xtensor-python: ${xtensor-python_INCLUDE_DIRS}")
endif()

find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

# Build
# =====

set(PY_FASTSCAPELIB _fastscapelib_py)

pybind11_add_module(${PY_FASTSCAPELIB}
    src/eroders.cpp
    src/flow_graph.cpp
    src/grid.cpp
    src/main.cpp)

target_link_libraries(${PY_FASTSCAPELIB}
                      PRIVATE fastscapelib xtensor-python Python::NumPy)

target_compile_features(${PY_FASTSCAPELIB} PUBLIC cxx_std_17)

target_compile_definitions(${PY_FASTSCAPELIB} PRIVATE VERSION_INFO=${PROJECT_VERSION})

# Installation
# ============

install(TARGETS ${PY_FASTSCAPELIB} LIBRARY DESTINATION .)
