// Copyright (c) 2024, twiinIT
//
// Proprietary license. All rights reserved.
//
// The full license is in the file LICENSE, distributed with this software.

#define FORCE_IMPORT_ARRAY

#include "block.hpp"
#include "flow_kernel.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pytensor.hpp"

#include <chrono>
#include <iostream>

namespace py = pybind11;

struct PyNumbaJitClass
{
    std::uintptr_t meminfoptr;
    std::uintptr_t dataptr;
};

struct PyNumbaFlowKernel
{
    std::uintptr_t func_ptr;
    std::uintptr_t node_data_getter_ptr;
    std::uintptr_t node_data_setter_ptr;
    std::uintptr_t node_data_create_ptr;
    std::uintptr_t node_data_free_ptr;
    PyNumbaJitClass data_ptr;
    int n_threads;
};


PYBIND11_MODULE(testpy, m)
{
    xt::import_numpy();

    py::class_<PyNumbaFlowKernel>(m, "Kernel")
        .def(py::init<>())
        .def_readwrite("func", &PyNumbaFlowKernel::func_ptr)
        .def_readwrite("node_data_getter", &PyNumbaFlowKernel::node_data_getter_ptr)
        .def_readwrite("node_data_setter", &PyNumbaFlowKernel::node_data_setter_ptr)
        .def_readwrite("node_data_create", &PyNumbaFlowKernel::node_data_create_ptr)
        .def_readwrite("node_data_free", &PyNumbaFlowKernel::node_data_free_ptr)
        .def_readwrite("data", &PyNumbaFlowKernel::data_ptr)
        .def_readwrite("n_threads", &PyNumbaFlowKernel::n_threads);

    py::class_<PyNumbaJitClass>(m, "JitClass")
        .def(py::init<>())
        .def_readwrite("meminfo", &PyNumbaJitClass::meminfoptr)
        .def_readwrite("data", &PyNumbaJitClass::dataptr);

    m.def("test_leak",
          [](PyNumbaFlowKernel py_kernel) -> int
          {
              // Release Python GIL
              py::gil_scoped_release release;
              // Cast pointers to opaque types
              auto kernel = (NumbaFlowKernel&) py_kernel;
              // Create a new node data
              NumbaJitClass new_node_data = kernel.node_data_create();

              return 0;
          });

    m.def("test_no_leak",
          [](PyNumbaFlowKernel py_kernel) -> int
          {
              // Release Python GIL
              py::gil_scoped_release release;
              // Cast pointers to opaque types
              auto kernel = (NumbaFlowKernel&) py_kernel;
              // Create a new node data
              NumbaJitClass new_node_data = kernel.node_data_create();
              // Then later delete it
              kernel.node_data_free(new_node_data.meminfoptr);

              return 0;
          });
}
