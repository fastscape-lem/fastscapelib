# Internals of Fastscapelib

This section provides information about Fastscapelib's internals, its design and
architecture. It is intended for anyone who wants to contribute improving
Fastscapelib or who wants to simply understand how it works under the hood.

:::{note}

Fastscapelib makes use of advanced techniques and patterns such as template
meta-programming, CRTP, type erasure, etc. Being familiar with these techniques
in C++ will help you better understand its internals.

:::

## Xtensor & Language Bindings

Fastscapelib heavily relies on the
[Xtensor](https://xtensor.readthedocs.io/en/latest/) C++ library for both its
API and its internals. Xtensor provides multi-dimensional array containers that
may be handled in a very similar way than [NumPy](https://numpy.org/) arrays in
Python.

For usage within C++ there are two main container types {cpp:type}`xt::xtensor`
and {cpp:type}`xt::xarray` with a static and dynamic number of dimensions,
respectively. In addition, Xtensor supports adapting other existing containers
and therefore enables seamless interoperability with array structures exposed in
other C++ libraries or other languages.

### Python

Fastscapelib's Python bindings are built with the help of
[Pybind11](https://pybind11.readthedocs.io) and
[Xtensor-python](https://xtensor-python.readthedocs.io/). The latter provides
two container types {cpp:type}`xt::pytensor` and {cpp:type}`xt::pyarray` that
each adapt CPython's [buffer
protocol](https://docs.python.org/3/c-api/buffer.html). This allows using the
Fastscapelib C++ algorithm implementations from Python directly with
[NumPy](https://numpy.org/) arrays without copying any data.

### Supporting Other Languages

Unlike for Python, 1st-class support for other language bindings is not planned
at the moment in Fastscapelib.

It shouldn't be too hard doing it, though, at least for a few other languages
popular in scientific computing. Like Xtensor-python for NumPy arrays,
[Xtensor-R](https://github.com/xtensor-stack/xtensor-r) and
[Xtensor-julia](https://github.com/xtensor-stack/xtensor-julia) provide
ready-to-use adaptors for using in-place R and Julia arrays, respectively.

Supporting languages like Fortran or MatLab would require writing [custom
Xtensor adaptors](https://xtensor.readthedocs.io/en/latest/adaptor.html).
Xtensor array containers support multiple layouts for representing
multi-dimensional arrays (i.e., row-major vs. column-major).

### Container Selectors

In order to avoid duplicating their implementation for each Xtensor container
type, most of Fastscapelib's classes (grid classes,
{cpp:class}`~fastscapelib::flow_graph`, etc.) expose a template parameter that
allows choosing the type of container to use for their array data members.

See {ref}`xtensor-selectors` and [Xtensor's
documentation](https://xtensor.readthedocs.io/en/latest/bindings.html#container-selection)
for more details.

## Grid Classes Hierarchy (CRTP)

All grid classes in Fastscapelib provide a common interface via their
{cpp:class}`~fastscapelib::grid` base class, which allows implementing algorithms
in a grid-agnostic fashion.

Grids may greatly differ from each other (e.g., 1-dimensional vs. 2-dimensional,
structured vs. unstructured) such that implementing a unified grid API is best
achieved using polymorphism.

Since the {cpp:class}`~fastscapelib::grid` API is often called a lot of times in
inner loops (i.e., iterate over grid nodes and their neighbors), static
polymorphism (at compile time) is preferred over dynamic polymorphism (at
runtime). It is implemented using the Curiously Recurring Template Pattern
([CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern)).

### Grid Inner Types

## Flow Graph

### Implementation(s)

### Python Bindings

Type erasure

## Flow Operators

Separate class vs. implementation.

Flow operator sequence (type erasure).
