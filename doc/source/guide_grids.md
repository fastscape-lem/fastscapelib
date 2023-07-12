(guide-grids)=
# Grids

Fastscapelib provides convenient classes for dealing with different kinds of
grids via some unified API, on top of which various algorithms (flow routing,
differential equation solvers, etc.) may be implemented in a grid-agnostic
fashion.

## Types of Grids

The table below describes the type of grids that are currently available in
Fastscapelib (also illustrated in {ref}`Figure 1 <fig_grid_types>`).

```{list-table} Grid types
:header-rows: 1
:name: table-grid-types
:widths: 25 20 55

* - C++ / Python Class
  - Properties
  - Usage Examples
* - {cpp:class}`~fastscapelib::profile_grid_xt` {py:class}`~fastscapelib.ProfileGrid`
  - 1-dimensional, uniform, static
  - Evolution of a river longitudinal profile, hillslope or escarpment cross-section.

    Solving flow routing or differential equations in this special case is often
    trivial, but this grid is still is useful for quickly switching between
    the 1D and 2D cases without the need to re-implement anything.
* - {cpp:class}`~fastscapelib::raster_grid_xt` {py:class}`~fastscapelib.RasterGrid`
  - 2-dimensional, rectangular, uniform, static
  - Evolution of an hillslope, escarpment, drainage basin, orogenic belt, etc.
* - {cpp:class}`~fastscapelib::trimesh_xt`
    {py:class}`~fastscapelib.TriMesh`
  - 2-dimensional, triangular (irregular) mesh, static
  - More complex than uniform grids but supports variable resolution
    (mesh refinement) and domain limits with complex shapes or planetary domains
    (no boundary).

    Single direction flow routing results in more "natural" river networks
    {cite:p}`Braun1997` than applied on raster grids (D8), although it is still
    constrained by the grid geometry.

    Note: this class doesn't implement any triangulation (triangles must be
    given to the grid constructor).
 ```

```{figure} _static/fig_grid_types.svg
:name: fig_grid_types
:alt: Grid types and representation
:align: center

Figure 1: Grid types (left: profile, center: raster, right: triangular mesh)
and representation (points: nodes, dashed lines: implicit or explicit
connectivity).
```

## Grid Representation

Fastscapelib currently implements a simple representation for all of its grids,
which only includes the grid nodes (points) and their connectivity ({ref}`Figure
1 <fig_grid_types>`).

Grid method names follow these conventions:

- ``nodes_xxx()``: get node-specific grid information either for all nodes, a
  subset of the nodes or one specific node depending on the (overloaded) method.
- ``neighbors_xxx()``: get grid information about the neighbors of a given grid
  node.

Although grid faces or cells are not represented explicitly (i.e., cell geometry
is not stored as grid data), some of their properties are exposed through the
grid API using the conventions above (e.g.,
{cpp:func}`~fastscapelib::grid::nodes_areas` returns the area of the direct
vicinity of each node).

Grid connectivity may be either implicit (e.g., in a
{cpp:class}`~fastscapelib::structured_grid` node neighbors are implied by the
grid structure) or explicit (e.g., the indices of triangle vertices must be
given to {cpp:type}`~fastscapelib::trimesh`). The
{cpp:func}`~fastscapelib::grid::neighbors` method and other related methods make
abstraction of the grid connectivity, which may be viewed as an implementation
detail.

See Sections {ref}`grid-node-iterators` and {ref}`connectivity-node-neighbors`
for more details on how to iterate through grid nodes and their neighbors.

(node-status-boundary-conditions)=
## Node Status and Boundary Conditions

Each node of the grid has a given status (or label), which usually serves to
determine how the model should behave at the domain limits, located on the grid
edges and/or inside the grid. All possible labels are defined in the
{cpp:enum}`~fastscapelib::node_status` (C++) and
{py:class}`~fastscapelib.NodeStatus` (Python) enum classes.

All grids expose a ``nodes_status`` parameter in their constructor, which
allows setting specific statuses for one or more nodes anywhere on the grid
(some restrictions may apply for the `LOOPED` status). Each grid also
exposes a ``nodes_status`` read-only getter or property that returns the
status of all grid nodes as an array.

Uniform grids {py:class}`~fastscapelib.ProfileGrid` and
{py:class}`~fastscapelib.RasterGrid` also provide convenient ways to set the
status of the edge or border nodes in their constructors, respectively via
{py:class}`~fastscapelib.ProfileBoundaryStatus` and
{py:class}`~fastscapelib.RasterBoundaryStatus` (which may be constructed
implicitly).

:::{note}

The status at nodes is defined once when creating a grid and cannot be changed
afterwards.

:::

:::{warning}

The usage or interpretation of grid node status and boundary conditions may
differ depending on the case. For example:

- {py:class}`~fastscapelib.FlowGraph` sets as {ref}`base levels
  <guide-base-level-nodes>` all nodes having the `FIXED_VALUE` status
  by default.

- {py:class}`~fastscapelib.DiffusionADIEroder` ignores the grid node status and
  always assumes fixed value boundaries on the raster grid borders.

:::

Below are a few examples of creating new grids with different statuses inside
the grid or on the grid boundaries. See also the {doc}`examples/index`.

1. A profile grid of 501 nodes with a total length equal to 500 m and with
   fixed-value boundaries, e.g., for simulating the evolution of a hillslope
   cross-section with base levels (river channel) at both sides.

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/grid/profile_grid.hpp"

namespace fs = fastscapelib;

fs::profile_boundary_status bs(fs::node_status::fixed_value);
auto grid = fs::profile_grid::from_length(501, 500.0, bs);
```

```{code-block} Python
:linenos:

import fastscapelib as fs

grid = fs.ProfileGrid.from_length(501, 500.0, NodeStatus.FIXED_VALUE)
```
````

2. A profile grid of 501 nodes with uniform node spacing of 1 meter, with
   ``CORE`` node status at both left and right edges and a fixed-value node at
   the middle of the grid, e.g., for simulating a valley cross-section:

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/grid/profile_grid.hpp"

namespace fs = fastscapelib;

fs::profile_boundary_status bs(fs::node_status::core);
auto grid = fs::profile_grid(501, 1.0, bs, { {250, fs::node_status::fixed_value} });
```

```{code-block} Python
:linenos:

import fastscapelib as fs

grid = fs.ProfileGrid(501, 1.0, fs.NodeStatus.CORE, {250: fs.NodeStatus.FIXED_VALUE})
```
````

3. A raster grid of 101x201 (row x col) nodes with a uniform node spacing of 100
   meters, with looped boundaries on the top-down borders, fixed-value on the
   left border and free (core) on the right border, e.g., for simulating an
   escarpment:

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

fs::raster_boundary_status bs{ fs::node_status::fixed_value,
                               fs::node_status::core,
                               fs::node_status::looped,
                               fs::node_status::looped };

auto grid = fs::raster_grid({ 101, 201 }, { 1e2, 1e2 }, bs);
```

```{code-block} Python
:linenos:

import fastscapelib as fs

bs = [
    fs.NodeStatus.FIXED_VALUE,
    fs.NodeStatus.CORE,
    fs.NodeStatus.LOOPED,
    fs.NodeStatus.LOOPED,
]

grid = fs.RasterGrid([101, 201], [1e2, 1e2], bs)
```
````

:::{tip}
How to choose the shape (number of nodes) of a uniform grid?

A round number + 1 like in the examples above results in a round value for the
uniform node spacing when the grid is constructed from a given, round total
length.

Alternatively, you might want to choose a power of two (e.g., 256, 512, 1024,
etc.) for optimal memory alignment of the grid data and internal arrays.

:::

## Grid Data Variables

Fastscapelib grid objects do not hold any grid data variable (fields), neither
by ownership nor by external reference (or pointer). This leaves users the
freedom of managing those data variables in a way that best suits their needs.

Grid data variables should be defined as
[Xtensor](https://xtensor.readthedocs.io)-compatible containers in C++ (e.g.,
{cpp:type}`xt::xtensor` or {cpp:type}`xt::xarray` objects) or NumPy arrays in Python.

The grid {cpp:func}`~fastscapelib::grid::size` and
{cpp:func}`~fastscapelib::grid::shape` methods are helpful for validating the
shape of the array containers at runtime.

Additionally, the {cpp:class}`~fastscapelib::grid_inner_types` C++ class
specialized for each grid exposes some properties like ``grid_data_type``,
``xt_selector`` and ``xt_ndims`` (for convenience, those are also exposed as
static data members of the grid classes). This is particularly useful for
maximizing code reusability in C++ (via type aliases) when creating new grid
data variables.

:::{note}

The dimension(s) of the grid data arrays (i.e., the length of
{cpp:func}`~fastscapelib::grid::shape` and the value of ``xt_dims``) do not
always match the number of physical axes of the grid. For example, variables on
a 2-dimensional triangular mesh are stored in 1-dimensional arrays.

:::

Below is an example of setting random initial topographic elevation values on a
raster grid.

````{tab-set-code}
```{code-block} C++
:linenos:

#include "xtensor/xarray.hpp"
#include "xtensor/xrandom.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

fs::raster_boundary_status bs{ fs::node_status::fixed_value };
fs::raster_grid grid({ 101, 101 }, { 200.0, 200.0 }, bs);

//
// setting a xt::xarray with the double data type
//
xt::xarray<double> elevation = xt::random::rand<double>(grid.shape());

//
// setting a xt::xtensor with dimensions and data type from grid inner types
//
using dtype = fs::raster_grid::grid_data_type;
using ndims = fs::grid_inner_types<fs::raster_grid>::xt_ndims;
xt::xtensor<dtype, ndims> elevation_alt = xt::random::rand<dtype>(grid.shape());

//
// use the xtensor container selector set for the grid
//
using selector = fs::raster_grid::xt_selector;
using ctype = fs::xt_container<selector, dtype, ndims>::tensor_type;
ctype elevation_alt2 = xt::random::rand<dtype>(grid.shape());
```

```{code-block} Python
:linenos:

import numpy as np
import fastscapelib as fs

grid = fs.RasterGrid([101, 101], [200.0, 200.0], fs.NodeStatus.FIXED_VALUE)

elevation = np.random.uniform(size=grid.shape)
```
````

(grid-node-iterators)=
## Grid Node Iterators

Fastscapelib provides a convenient, stl-compatible way to iterate over the grid
nodes in C++, possibly filtered by {ref}`node status
<node-status-boundary-conditions>`, using
{cpp:func}`~fastscapelib::grid::nodes_indices`.

:::{note}

Iterating over grid nodes using pure-Python loops is very slow. In general it is
possible to obtain the same results much faster using operations on NumPy arrays
(vectorization). If that is not possible, great speed-ups can still be obtained,
e.g., by using [Numba](http://numba.pydata.org/)'s jit compiler (non-object
mode).

:::

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

using fixed_value = fs::node_status::fixed_value;
fs::raster_boundary_status bs(fixed_value);
fs::raster_grid grid({ 101, 101 }, { 200.0, 200.0 }, bs);

//
// iterate over all grid nodes
//
for (auto& idx_flat : grid.nodes_indices())
{
    // on a raster grid, flat index corresponds to the
    // current row index * total nb. of rows + current column index
    std::cout << "current node flat index: " << idx_flat << std::endl;
}

//
// iterate over fixed-value nodes (all border nodes in this case)
//
for (auto& idx_flat : grid.nodes_indices(fixed_value))
{
    std::cout << "current node flat index: " << idx_flat << std::endl;
}

```
```{code-block} Python
:linenos:

import fastscapelib as fs
import numpy as np

fixed_value = fs.NodeStatus.FIXED_VALUE
grid = fs.RasterGrid([100, 100], [200.0, 200.0], fixed_value)

# iterate over all grid nodes (slow!)
for idx_flat in grid.nodes_indices():
    print(f"current node flat index: {idx_flat}")

# iterate over fixed-value nodes (slow!)
for idx_flat in grid.nodes_indices(fixed_value):
    print(f"current node flat index: {idx_flat}")

# do vectorized operations when possible (fast!)
data = np.random.rand(size=grid.shape)
data.ravel()[grid.node_indices(fixed_value)] = 0.0

```
````

(connectivity-node-neighbors)=
## Connectivity and Node Neighbors

Iterating over the neighbors of a grid node is simple, does not depend on the
type of grid and follows looped boundaries if any.

A typical way to iterate over grid nodes and their neighbors is illustrated in
the example below.

:::{note}

Iterating over grid nodes and their neighbors in Python is slow. Grid neighbors
methods are still exposed in Python for development and debugging purposes.

Note that it is not possible to use Fastscapelib grid (Python) objects directly
in [Numba](http://numba.pydata.org/) jitted functions in non-object mode.

:::


````{tab-set-code}
```{code-block} C++
:linenos:

#include <iostream>
#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

fs::raster_boundary_status bs(fs::node_status::fixed_value);
fs::raster_grid grid({ 101, 101 }, { 200.0, 200.0 }, bs);

//
// optimization: initialize the neighbors container out of the loops
//
fs::raster_grid::neighbors_type neighbors;

for (auto& idx_flat : grid.nodes_indices())
{
    for (auto& nb : grid.neighbors(idx_flat, neighbors))
    {
        std::cout << "flat index: " << nb.idx << std::endl;
        std::cout << "distance to neighbor: " << nb.distance << std::endl;
        std::cout << "status of neighbor: " << nb.status << std::endl;
    }
}
```
```{code-block} Python
:linenos:

import fastscapelib as fs

grid = fs.RasterGrid([100, 100], [200.0, 200.0], fs.NodeStatus.FIXED_VALUE)

for idx_flat in range(grid.size):
    neighbors = grid.neighbors(idx_flat)
    for nb in neighbors:
        print(f"flat index: {nb.idx}")
        print(f"distance to neighbor: {nb.distance}")
        print(f"status of neighbor: {nb.status}")
```
````

The example above is using the simple {cpp:class}`~fastscapelib::neighbor` (C++)
and {py:class}`~fastscapelib.Neighbor` (Python) structures to store information
about a node neighbor. If you just need the total number of neighbors for a
given grid node or only the neighbor (flat) indices or the distances from one
node to its neighbors, you can use the alternative methods
``neighbors_count()``, ``neighbors_indices()`` and ``neighbors_distances()``.

{cpp:class}`~fastscapelib::raster_grid_xt` (C++) also provides some method
overloads to use row and column indices instead of grid flat indices.

### Raster Grid Connectivity

{cpp:class}`~fastscapelib::raster_grid_xt` (C++) exposes a template parameter
that allows choosing the grid connectivity among three modes (see
{cpp:enum}`~fastscapelib::raster_connect` and {ref}`Figure 2
<fig_raster_connect>`):

- ``queen`` (8-node connectivity, including diagonals, set by default)
- ``rook`` (4-node connectivity, no diagonal)
- ``bishop`` (4-nodes connectivity, diagonals only)

:::{note}

{py:class}`~fastscapelib.RasterGrid` (Python) only supports the ``queen`` mode.

:::

```{figure} _static/fig_raster_connect.svg
:name: fig_raster_connect
:alt: Raster connectivity modes
:align: center

Figure 2: Raster grid connectivity modes: queen (left), rook (center) and bishop (right).
```

### Caching Neighbor Node Indices

Fastscapelib implements a cache system for retrieving neighbor node indices,
which is particularly useful for speeding-up neighbor lookup on structured
(raster) grids. Such a cache helps avoiding repetitive operations like checking
if the current visited node is on the boundary or finding neighbors on looped
boundaries. Caveat: this speed-up is achieved at the expense of memory
consumption.

The cache is exposed as a template parameter for
{cpp:class}`~fastscapelib::raster_grid_xt` so that this behavior may be
customized (cache may be disabled). For other grid types like
{cpp:class}`~fastscapelib::trimesh_xt` with explicit topology, the
cache is not needed and not used.

See {cpp:class}`~fastscapelib::neighbors_cache` and
{cpp:class}`~fastscapelib::neighbors_no_cache` for more details.

:::{note}

In the Python bindings the cache is enabled for raster grids and there is
currently no option to disable it.

:::
