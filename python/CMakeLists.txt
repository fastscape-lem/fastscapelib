cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH})

# Dependencies
# ============

# -- Build Python bindings independently of fastscapelib
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(fastscapelib-python)

    find_package(xtensor REQUIRED CONFIG)
    find_package(fastscapelib REQUIRED CONFIG)
endif()

# -- Other deps
find_package(xtensor-python REQUIRED)
message(STATUS "Found xtensor-python: ${xtensor-python_INCLUDE_DIRS}/xtensor-python")

find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

find_package(NumPy REQUIRED)
message(STATUS "Found NumPy: v${NUMPY_VERSION} ${NUMPY_INCLUDE_DIRS}")

# Installation directory
# ======================

# -- Set the Python site-package installation directory
find_package(PythonInterp)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; import os;print(get_python_lib())"
                OUTPUT_VARIABLE PY_SP_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

set(SP_DIR "${PY_SP_DIR}" CACHE PATH "Site-package directory")
message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> site-packages: " ${SP_DIR})

# Build
# =====

set(PY_FASTSCAPELIB _fastscapelib_py)

pybind11_add_module(${PY_FASTSCAPELIB}
    src/algo.cpp
    src/eroders.cpp
    src/flow_graph.cpp
    src/grid.cpp
    src/main.cpp)

target_include_directories(${PY_FASTSCAPELIB} PUBLIC ${NUMPY_INCLUDE_DIRS})

target_link_libraries(${PY_FASTSCAPELIB}
                      PUBLIC fastscapelib xtensor-python)

target_compile_features(${PY_FASTSCAPELIB} PUBLIC cxx_std_17)

# Installation
# ============

install(TARGETS ${PY_FASTSCAPELIB}
    DESTINATION ${SP_DIR})

# install type annotations
install(FILES _fastscapelib_py.pyi
    DESTINATION ${SP_DIR})
