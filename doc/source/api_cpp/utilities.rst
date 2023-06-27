Utilities
=========

.. _xtensor-selectors:

Xtensor container selectors
---------------------------

Defined in ``fastscapelib/utils/xtensor_utils.hpp``

Most Fastscapelib data structures (i.e., grids, flow graph, etc.) expose a
template parameter that allows setting the type of xtensor container to use for
their public array members.

This allows reusing the same classes (and avoid conversion copies) in various
contexts, e.g.,

- in C++ applications with :cpp:type:`xt::xtensor` / :cpp:type:`xt::xarray`
- in Python bindings with :cpp:type:`xt::pytensor` / :cpp:type:`xt::pyarray`
- etc.

|

.. doxygenstruct:: fastscapelib::xt_selector

.. doxygenstruct:: fastscapelib::xt_container

.. doxygentypedef:: fastscapelib::xt_tensor_t

.. doxygentypedef:: fastscapelib::xt_array_t

|


Iterators and virtual containers
--------------------------------

Defined in ``fastscapelib/utils/iterators.hpp``

For convenience, Fastscapelib provides STL-compatible iterators and virtual
containers for looping over grid or flow graph nodes.

|

.. doxygenclass:: fastscapelib::grid_nodes_indices

.. doxygenclass:: fastscapelib::stl_container_iterator_wrapper
