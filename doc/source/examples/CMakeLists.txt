cmake_minimum_required(VERSION 3.28)

project(test)

set(CMAKE_CXX_STANDARD 20)
add_executable(test main.cpp)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

find_package(pybind11 CONFIG REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-python REQUIRED)

pybind11_add_module(testpy
    main.cpp
)
target_link_libraries(testpy PRIVATE xtensor xtensor-python Python::NumPy)
