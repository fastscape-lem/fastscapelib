(guide-grids)=
# Grids

Fastscapelib provides convenient classes for dealing with different kinds of
grids via some unified API, on top of which various algorithms (flow routing,
differential equation solvers, etc.) may be implemented in a grid-agnostic
fashion.

## Types of Grids

The table below describes the type of grids that are currently available in Fastscapelib.

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

    Algorithm implementations on this grid are often trivial, but it is useful for
    quickly switching between the 1D and 2D cases while experimenting.
* - {cpp:class}`~fastscapelib::raster_grid_xt` {py:class}`~fastscapelib.RasterGrid`
  - 2-dimensional, rectangular, uniform, static
  - Evolution of an hillslope, escarpment, drainage basin, orogenic belt, etc.
* - {cpp:class}`~fastscapelib::unstructured_mesh_xt`
    {py:class}`~fastscapelib.UnstructuredMesh`
  - 2-dimensional, irregular (triangulated) mesh, static
  - More complex than uniform grids but supports variable resolution
    (mesh refinement) and domain limits with complex shapes or planetary domains
    (no boundary).

    Single direction flow routing results in more "natural" river networks
    {cite:p}`Braun1997` than applied on raster grids (D8).

    Note: this class doesn't implement any triangulation (triangles must be
    provided).
 ```

## Grid Representation

Fastscapelib currently implements a simple representation for all of its grids,
which only includes the grid nodes (points) and their implicit or explicit
connectivity. See Sections {ref}`grid-node-iterators` and
{ref}`connectivity-node-neighbors` for how to iterate through grid nodes and
their neighbors.

Grid faces or cells are not represented explicitly, i.e., cell vertices and
edges are not stored anywhere. However, the area of the cell surrounding a given
node can be accessed via the {cpp:func}`~fastscapelib::grid::node_area` method
(C++ only).

(node-status-boundary-conditions)=
## Node Status and Boundary Conditions

Each node of the grid has a given status (or label), which usually serves to
determine how the model should behave at the domain limits, located on the grid
edges and/or inside the grid. All possible labels are defined in the
{cpp:enum}`~fastscapelib::node_status` (C++) and
{py:class}`~fastscapelib.NodeStatus` (Python) enum classes.

All grids expose a ``status_at_nodes`` parameter in their constructor, which
allows setting specific statuses for one or more nodes anywhere on the grid
(some restrictions may apply for the `LOOPED_BOUNDARY` status). Each grid also
exposes a ``status_at_nodes`` read-only property that returns the status of all
grid nodes as an array.

Uniform grids {py:class}`~fastscapelib.ProfileGrid` and
{py:class}`~fastscapelib.RasterGrid` also provide convenient ways to set the
status of the edge or border nodes in their constructors, respectively using
{py:class}`~fastscapelib.ProfileBoundaryStatus` and
{py:class}`~fastscapelib.RasterBoundaryStatus`.

:::{note}

The status at nodes is defined once when creating a grid and cannot be changed
afterwards.

:::

:::{warning}

The usage or interpretation of grid node status and boundary conditions may
differ from one component to another. For example:

- {py:class}`~fastscapelib.FlowGraph` sets as {ref}`base levels
  <guide-base-level-nodes>` all nodes having the `FIXED_VALUE_BOUNDARY` status
  by default.

- {py:class}`~fastscapelib.DiffusionADIEroder` ignores the grid node status and
  always assume fixed value boundaries on the raster grid borders.

:::

Below are a few examples of creating new grids with different statuses inside
the grid or on the grid boundaries. See also the {doc}`examples/index`.

1. A profile grid of 501 nodes (including edge nodes), a total length equal to
   500 m and with fixed-value boundaries, e.g., for simulating the evolution of
   a hillslope cross-section with base levels (river channel) at both sides.

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/grid/profile_grid.hpp"

namespace fs = fastscapelib;

//
// fs::profile_boundary_status is constructed implicitly from the 3rd parameter.
//
auto grid = fs::profile_grid::from_length(
    501, 500.0, fs::node_status::fixed_value_boundary);
```

```{code-block} Python
:linenos:

import fastscapelib as fs

boundaries = fs.ProfileBoundaryStatus(fs.NodeStatus.FIXED_VALUE_BOUNDARY)
grid = fs.ProfileGrid.from_length(501, 500.0, boundaries, [])
```
````

2. A profile grid of 501 nodes (including edge nodes), a uniform node spacing of
   1 meter, with ``CORE`` node status at both left and right edges and a
   fixed-value node in the middle, e.g., for simulating a valley cross-section
   with a single river (base level):

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/grid/profile_grid.hpp"

namespace fs = fastscapelib;

auto grid = fs::profile_grid(
    501,
    1.0,
    fs::node_status::core,
    { fs::node({ 250, fs::node_status::fixed_value_boundary }) });
```

```{code-block} Python
:linenos:

import fastscapelib as fs

boundaries = fs.ProfileBoundaryStatus(fs.NodeStatus.CORE)
grid = fs.ProfileGrid(
    501,
    500.0,
    boundaries,
    [fs.Node(250, fs.NodeStatus.FIXED_VALUE_BOUNDARY)],
)
```
````

3. A raster grid of 101x201 (row x col) nodes, a uniform node spacing of 100
   meters, with looped boundaries on the top-down borders, fixed-value on the
   left border and free (core) on the right border, e.g., for simulating an
   escarpment:

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

fs::raster_boundary_status bs{ fs::node_status::fixed_value_boundary,
                               fs::node_status::core,
                               fs::node_status::looped_boundary,
                               fs::node_status::looped_boundary };

auto grid = fs::raster_grid({ 101, 201 }, { 1e2, 1e2 }, bs);
```

```{code-block} Python
:linenos:

import fastscapelib as fs

bs = fs.RasterBoundaryStatus(
    [
        fs.NodeStatus.FIXED_VALUE_BOUNDARY,
        fs.NodeStatus.CORE,
        fs.NodeStatus.LOOPED_BOUNDARY,
        fs.NodeStatus.LOOPED_BOUNDARY,
    ]
)

grid = fs.RasterGrid([101, 201], [1e2, 1e2], bs, [])
```
````

## Grid Data Variables

Fastscapelib's grid objects do not hold any grid data variable (fields), neither
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
a 2-dimensional unstructured mesh are stored in 1-dimensional arrays.

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

fs::raster_boundary_status boundaries{ fs::node_status::fixed_value_boundary };
fs::raster_grid grid({ 101, 101 }, { 200.0, 200.0 }, boundaries);

//
// setting a xt::xarray with double data type
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

boundaries = fs.RasterBoundaryStatus(fs.NodeStatus.FIXED_VALUE_BOUNDARY)
grid = fs.RasterGrid([101, 101], [200.0, 200.0], boundaries, [])

elevation = np.random.uniform(size=grid.shape)
```
````

(grid-node-iterators)=
## Grid Node Iterators

Fastscapelib provides a convenient, stl-compatible way to iterate over grid
nodes in C++, possibly filtered by {ref}`node status
<node-status-boundary-conditions>`, using
{cpp:func}`~fastscapelib::grid::node_indices`.

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

using fixed_value = fs::node_status::fixed_value_boundary;
fs::raster_grid grid({ 101, 101 }, { 200.0, 200.0 }, fixed_value);

//
// iterate over all grid nodes
//
for (auto& idx_flat : grid.node_indices())
{
    // on a raster grid, flat index corresponds to the
    // current row index * total nb. of rows + current column index
    std::cout << "current node flat index: " << idx_flat << std::endl;
}

//
// iterate over fixed-value nodes (all border nodes in this case)
//
for (auto& idx_flat : grid.node_indices(fixed_value))
{
    std::cout << "current node flat index: " << idx_flat << std::endl;
}

```
```{code-block} Python
:linenos:

import fastscapelib as fs

boundaries = fs.RasterBoundaryStatus(fs.NodeStatus.FIXED_VALUE_BOUNDARY)
grid = fs.RasterGrid([100, 100], [200.0, 200.0], boundaries, [])

node_status_flat = grid.status_at_nodes.ravel()

# iterate over all grid nodes
for idx_flat in range(grid.size):
    print(f"current node flat index: {idx_flat}")

# iterate over fixed-value nodes (all border nodes in this case)
for idx_flat in range(grid.size):
    if node_status_flat[idx_flat] == fs.NodeStatus.FIXED_VALUE_BOUNDARY:
        print(f"current node flat index: {idx_flat}")
```
````

(connectivity-node-neighbors)=
## Connectivity and Node Neighbors

Iterating over the neighbors of a node is simple, does not depend on the type
of grid and/or boundary conditions, and follows looped boundaries if any.

A typical way to iterate over grid nodes and their neighbors is as follows:

````{tab-set-code}
```{code-block} C++
:linenos:

#include <iostream>
#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

using fixed_value = fs::node_status::fixed_value_boundary;
fs::raster_grid grid({ 101, 101 }, { 200.0, 200.0 }, fixed_value);

//
// initialize the neighbors container out of the iteration loops
// for efficiency
//
fs::raster_grid::neighbors_type neighbors;

for (auto& idx_flat : grid.node_indices())
{
    for (auto& nb : grid.neighbors(idx, neighbors))
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

boundaries = fs.RasterBoundaryStatus(fs.NodeStatus.FIXED_VALUE_BOUNDARY)
grid = fs.RasterGrid([100, 100], [200.0, 200.0], boundaries, [])

for idx_flat in range(grid.size):
    neighbors = grid.neighbors(idx_flat)
    for nb in neighbors:
        print(f"flat index: {nb.idx}")
        print(f"distance to neighbor: {nb.distance}")
        print(f"status of neighbor: {nb.status}")
```
````

:::{note}

Iterating over grid nodes and their neighbors in Python is slow. Grid neighbors
methods are still exposed in Python for development and debugging purposes.

:::

The example above is using the simple {cpp:class}`~fastscapelib::neighbor` (C++)
and {py:class}`~fastscapelib.Neighbor` (Python) structures to store information
about a node neighbor. Sometimes we just need the total number of neighbors for
a given node or only the neighbor (flat) indices or the distances from one node
to its neighbors. There are alternative methods for that, i.e.,
``neighbors_count()``, ``neighbors_indices()`` and ``neighbors_distances()``.

{cpp:class}`~fastscapelib::raster_grid_xt` also provides some method overloads
to use row and column indices instead of flat grid indices.

### Raster Grid Connectivity

### Caching Neighbor Node Indices

Fastscapelib implements a caching system for retrieving neighbor node indices,
which is particularly useful for speeding-up neighbor lookup on structured
(raster) grids. Such a cache helps avoiding repetitive operations like checking
if the current is on the boundary and finding reflective neighbors (looped
boundary) if any. However, that speed-up is achieved at the expense of memory.

The cache is exposed as a template parameter for
{cpp:class}`~fastscapelib::raster_grid_xt` so that this behavior may be
customized (cache may be disabled). For other grid types like
{cpp:class}`~fastscapelib::unstructured_mesh_xt` with explicit topology, the
cache is not needed and not used.

See {cpp:class}`~fastscapelib::neighbors_cache` and
{cpp:class}`~fastscapelib::neighbors_no_cache` for more details.

:::{note}

In the Python bindings the cache is enabled for raster grids and there is
currently no option to disable it.

:::
