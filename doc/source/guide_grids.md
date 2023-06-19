(guide-grids)=
# Grids

## Types of Grids

## Node Status and Boundary Conditions

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

Additionally, the {cpp:class}`~fastscapelib::grid_inner_types` structure
specialized for each grid exposes some properties like ``grid_data_type``,
``xt_selector`` and ``xt_ndims`` (for convenience, those are also exposed as
static data members of the grid classes). This is particularly useful for
maximizing code reusability in C++ (via type aliases) when creating new grid
data variables.

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
fs::raster_grid grid({ 100, 100 }, { 200.0, 200.0 }, boundaries);

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
grid = fs.RasterGrid([100, 100], [200.0, 200.0], boundaries, [])

elevation = np.random.uniform(size=grid.shape)
```
````

## Grid Node Iterators

## Connectivity and Node Neighbors

:::{note}

Iterating over grid nodes and their neighbors in Python is slow. Grid neighbors
methods are still exposed in Python for development and debugging purpose.

:::

### Caching

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

:::{note}

In the Python bindings the cache is enabled for raster grids and there is
currently no option to disable it

:::

See {cpp:class}`~fastscapelib::neighbors_cache` and
{cpp:class}`~fastscapelib::neighbors_no_cache` for more details.
