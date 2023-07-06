(internals)=
# Internals of Fastscapelib

This section provides information about Fastscapelib's internals, its design and
architecture. It is intended for anyone who wants to contribute improving
Fastscapelib or who wants to simply understand how it works under the hood.

:::{note}

Fastscapelib makes use of advanced techniques and patterns such as template
meta-programming, CRTP, type erasure, etc. Being familiar with these techniques
in C++ will help you better understand its internals.

:::

(internals-xtensor-language-bindings)=
## Xtensor & Language Bindings

Fastscapelib heavily relies on the [Xtensor](https://xtensor.readthedocs.io) C++
library for both its API and its internals. Xtensor provides multi-dimensional
array containers that may be handled in a very similar way than
[NumPy](https://numpy.org/) arrays in Python.

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

Unlike for Python, first-class support for other language bindings is not
planned at the moment.

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

## Grids

### Class Hierarchy (CRTP)

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
([CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern)):

- the {cpp:class}`~fastscapelib::grid` and
  {cpp:class}`~fastscapelib::structured_grid` base classes both have a ``G``
  template parameter that represents one of the grid final (leaf) classes
- those base classes also have a ``derived_grid()`` protected method to cast the
  grid instance into the ``G`` grid (leaf) type

See the class diagram below (showing only a subset of the public API for
clarity). See also the {doc}`grid API reference <api_cpp/grid>` for more
details.

```{mermaid}
classDiagram
    grid~G~ <|-- structured_grid~G~
    grid~G~ <|-- trimesh_xt
    structured_grid~G~ <|-- profile_grid_xt
    structured_grid~G~ <|-- raster_grid_xt
    class grid~G~{
        #derived_grid()
        +is_structured()
        +is_uniform()
        +neighbors_count(idx)
        +neighbors_indices(idx)
        +neighbors_distances(idx)
        +neighbors(idx)
    }
    class structured_grid~G~{
        +spacing()
        +length()
        +shape()
    }
    class profile_grid_xt
    class raster_grid_xt{
        +neighbors_indices(row, col)
        +neighbors(row, col)
    }
    class trimesh_xt
```

### Python Bindings

Note that since we use static polymorphism, only the grid leaf classes are
exposed in Python, i.e., {cpp:class}`~fastscapelib::grid` and
{cpp:class}`~fastscapelib::structured_grid` have no Python bindings and the grid
Python classes have no common base class.

### Grid Inner Types

The grid base classes need some types and values that are specific to each grid
derived type (e.g., grid resolution, number of dimensions of grid field arrays,
maximum number of node neighbors, etc.) and that must therefore be defined in
grid leaf classes.

However, with CRTP this is not possible since a CRTP leaf class is only declared
when the CRTP base class is being defined.

As a workaround, those grid inner types and values are defined in a separate
template class {cpp:class}`~fastscapelib::grid_inner_types` that is specialized
for each grid leaf class.

## Flow Graph

(internals-flow-graph-impl)=

### Implementation(s)

The template class {cpp:class}`~fastscapelib::flow_graph` is distinct from its
implementation done in the template class `flow_graph_impl<G, S, Tag>`. Both
expose the same template parameters, where ``Tag`` is used for selecting a given
graph implementation (or representation).

Currently, Fastscapelib only supports one implementation based on fixed-shape
arrays to store the graph topology. This implementation is selected via the
{cpp:class}`~fastscapelib::flow_graph_fixed_array_tag` type. Other
implementations may be added in the future, for example based on linked-lists or
sparse matrices, maybe reusing C++ graph libraries like [the Boost Graph
Library](https://www.boost.org/doc/libs/1_82_0/libs/graph/doc/index.html) or a
[GraphBLAS](https://graphblas.org/) C++ API.

The same template parameter ``Tag`` is also used by the {ref}`flow operator
implementation classes <internals-flow-operators>`.

### Python Bindings

It is not possible with [Pybind11](https://pybind11.readthedocs.io) to expose in
Python a C++ template class that is not fully specialized.

Since the template class {cpp:class}`~fastscapelib::flow_graph` has the grid
(leaf) type as template parameter ``G``, we rely on the [type
erasure](https://en.wikipedia.org/wiki/Type_erasure) technique in order to avoid
exposing separate flow graph classes in Python for each grid type.

(internals-flow-operators)=

## Flow Operators

Like the flow graph, the architecture of flow operators decouples an operator
from its implementation(s) so that it will be possible to later support other
graph internal representations while reusing the same operators. It also makes
easier exposing the flow operator classes in Python, since
{cpp:class}`~fastscapelib::flow_operator` is not a template class, nor any of
its subclasses.

The architecture of flow operators is a bit more complex, as shown in the
diagram below.

```{mermaid}
classDiagram
    flow_operator "1" *-- "1" flow_operator_impl_base~FG, OP~: shared pointer
    flow_operator_impl_base~OP~ <|-- flow_operator_impl~FG, OP, Tag~
    flow_operator_impl~FG, OP, Tag~ "1" *-- "1" flow_operator_impl_facade~FG~ : type erasure
    flow_operator_sequence~FG~ "1" --* "n" flow_operator_impl_facade~FG~
    flow_operator_sequence~FG~ "1" *-- "1" flow_graph~G, S, Tag~
    class flow_operator{
        +elevation_updated
        +graph_updated
        +in_flowdir
        +out_flowdir
        +virtual name()
    }
    class flow_graph~G, S, Tag~{
        +operators()
    }
    class flow_operator_impl_base~FG, OP~{
        #m_op_ptr
        +apply(...)
        +save(...)
    }
    class flow_operator_impl~FG, OP, Tag~
    class flow_operator_impl_facade~FG~{
        +apply(...)
        +save(...)
    }
    class flow_operator_sequence~FG~{
        +elevation_updated()
        +graph_updated()
        +out_flowdir()
    }
```

### Operator (Interface) Classes

The {cpp:class}`~fastscapelib::flow_operator` abstract base class provides the
operator interface. Each operator must be defined as a subclass, within which
the operator parameters (if any) should also be defined.

### Operator Implementation (Template) Class

The ``flow_operator_impl<FG, OP, Tag>`` template class is where all the
operators are implemented. It must be partially specialized for the two
following template parameters:

- ``OP`` (flow operator type): at least one specialization must exist for each
  {cpp:class}`~fastscapelib::flow_operator` subclass.
- ``Tag`` ({ref}`flow graph implementation tag <internals-flow-graph-impl>`): a
  flow operator may be implemented for one or more of the available flow graph
  representations.

It is not specialized regarding the ``FG`` parameter (actual flow graph
implementation type), which also depends on the grid (leaf) type ``G`` passed
from the flow graph (i.e., ``FG`` is equivalent to ``flow_graph_impl<G, S,
Tag>``).

``flow_operator_impl<FG, OP, Tag>`` inherits from a base class that has the
following methods and members:

- ``apply()``: method through which the actual logic of the operator is executed
  and that should be re-implemented by every operator
- ``save()``: method re-implemented by {cpp:class}`~fastscapelib::flow_snapshot`
  and generally not re-implemented by other operators
- ``m_op_ptr``: a shared pointer to an instance of ``OP``, useful for accessing
  from within the implementation code the current values set for the operator
  parameters, if any

### Operator Sequence Class

The template class ``flow_operator_sequence<FG>`` serves as an intermediate for
instantiating the implementation of each of the operators that have been passed
to {cpp:class}`~fastscapelib::flow_graph`. This is done via the template class
``flow_operator_impl_facade<FG>``, which implements [type
erasure](https://en.wikipedia.org/wiki/Type_erasure) so that instances of
``flow_operator_impl<FG, OP, Tag>`` of arbitrary operators types ``OP`` can be
built and added to the sequence at runtime.
