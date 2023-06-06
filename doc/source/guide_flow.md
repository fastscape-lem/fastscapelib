# Flow Routing

Flow routing consists of finding the paths along which some material (e.g.,
water, sediment) "flows" over the topographic surface, from upslope to downslope
areas. It is an important component of most Landscape Evolution Models (LEMs) ;
some of the main driver processes like river channel erosion depend on it.

Flow routing is also a complex problem that requires advanced algorithms and
(graph-based) data structures in order to compute the flow paths more-or-less
accurately, which often represents the main bottleneck in landscape evolution
simulations.

Fastscapelib implements a collection of efficient, state-of-the-art algorithms
for flow routing (e.g., {cite:p}`Braun2013`, {cite:p}`Barnes2014`,
{cite:p}`Cordonnier2019`). Those algorithms are optimized for LEMs (i.e.,
require updating the flow paths repeatedly) and may be run through a flexible
and user-friendly API.

## Flow Graph

Flow paths are represented explicitly in a "flow graph" data structure (more
precisely, a Directed Acyclic Graph or DAG) that connects each grid node to its
flow receiver(s) (neighbors).

Although it can be handled as a standalone structure, a flow graph is closely
related to a grid (Figure TODO):

- there is a perfect match between the sets of graph nodes and grid nodes
- the set of graph edges corresponds to a subset of the grid (implicit or
  explicit) edges

A {cpp:class}`~fastscapelib::flow_graph` (C++) or
{py:class}`~fastscapelib.FlowGraph` (Python) object can be created from any grid
object (see Section TODO), e.g., in the example below from a raster grid:

````{tab-set-code}
```{code-block} C++
:linenos:
:emphasize-lines: 10

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

fs::raster_boundary_status boundaries{ fs::node_status::fixed_value_boundary };
fs::raster_grid grid({ 100, 100 }, { 200.0, 200.0 }, boundaries);

fs::flow_graph<fs::raster_grid> graph(grid, { fs::single_flow_router() });
```

```{code-block} Python
:linenos:
:emphasize-lines: 6

import fastscapelib as fs

boundaries = fs.RasterBoundaryStatus(fs.NodeStatus.FIXED_VALUE_BOUNDARY)
grid = fs.RasterGrid([100, 100], [200.0, 200.0], boundaries, [])

graph = fs.FlowGraph(grid, [fs.SingleFlowRouter()])
```
````

A sequence of flow operators must be passed as second argument to the flow graph
constructor in order to define the flow routing strategy (see Section TODO
below).

### Base Level Nodes

The flow routing algorithms implemented in Fastscapelib require one or more
nodes of the graph to be set as "base level" (or outlet) nodes.

These specific nodes may collect input flow but cannot release any output flow
(or that have an unknown flow receiver).

By default, all grid nodes with status
{py:attr}`~fastscapelib.NodeStatus.FIXED_VALUE_BOUNDARY` are set as base
level nodes when creating a flow graph.

### Initializing / Updating the Flow Paths

Creating a new flow graph object doesn't compute any flow path yet. To
(re)compute the flow paths, the {py:meth}`~fastscapelib.FlowGraph.update_routes`
method must be called with an input topographic surface (i.e., an elevation
field defined on the grid):

````{tab-set-code}
```{code-block} C++
:linenos:
:emphasize-lines: 17

#include "xtensor/xarray.hpp"
#include "xtensor/xrandom.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

fs::raster_boundary_status boundaries{ fs::node_status::fixed_value_boundary };
fs::raster_grid grid({ 100, 100 }, { 200.0, 200.0 }, boundaries);

fs::flow_graph<fs::raster_grid> graph(grid, { fs::single_flow_router() });

xt::xarray<double> elevation = xt::random::rand<double>(grid.shape());

// update_routes returns a const reference!
const auto new_elevation = graph.update_routes(elevation);
```

```{code-block} Python
:linenos:
:emphasize-lines: 11

import numpy as np
import fastscapelib as fs

boundaries = fs.RasterBoundaryStatus(fs.NodeStatus.FIXED_VALUE_BOUNDARY)
grid = fs.RasterGrid([100, 100], [200.0, 200.0], boundaries, [])

graph = fs.FlowGraph(grid, [fs.SingleFlowRouter()])

elevation = np.random.uniform(size=grid.shape)

new_elevation = graph.update_routes(elevation)
```
````

{py:meth}`~fastscapelib.FlowGraph.update_routes` returns another elevation
field, which may or may not differ from the input elevation field depending on
the applied flow operators. Some operators fill the closed depressions found in
the input topography.

:::{note}

{py:meth}`~fastscapelib.FlowGraph.update_routes` never updates in-place the
values of the elevation field passed as input.

 :::

## Flow Operators

Flow can be routed over the topographic surface in many different ways ;
choosing one approach over another highly depends on the case studied.
Fastscapelib relies on the concept of "flow operators" that provide a flexible
and convenient solution for implementing simple to advanced flow routing
strategies.

A flow operator is a "routing unit" that:

- may read and/or modify in-place the flow graph instance to which it has been attached.
- may read and/or update (in a copy) the elevation values passed to the
  {py:meth}`~fastscapelib.FlowGraph.update_routes` method.
- may expose zero, one or more parameters with values that can be changed at any
  time between two consecutive calls to
  {py:meth}`~fastscapelib.FlowGraph.update_routes`.

There are currently three categories of operators (see the {ref}`C++
<cpp-api-flow-operators>` and {ref}`Python <py-api-flow-operators>` API
reference for a full list of available operators).

### Flow Routers

Flow router operators generally compute new flow paths from scratch and fully
update the flow graph. The updated flow graph has a defined
{py:class}`~fastscapelib.FlowDirection`: either ``SINGLE`` (one unique receiver
per node) or ``MULTI`` (flow partitioned among multiple receiver nodes with
fixed or variable weights).


````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/flow/flow_router.hpp"

namespace fs = fastscapelib;

fs::single_flow_router::graph_updated  // true
fs::single_flow_router::out_flowdir    // fs::flow_direction::single

fs::multi_flow_router::graph_updated   // true
fs::multi_flow_router::out_flowdir     // fs::flow_direction::multi
```

```{code-block} Python
:linenos:

import fastscapelib as fs

fs.SingleFlowRouter.graph_updated   # True
fs.SingleFlowRouter.out_flowdir     # fs.FlowDirection.SINGLE

fs.MultiFlowRouter.graph_updated    # True
fs.MultiFlowRouter.out_flowdir      # fs.FlowDirection.MULTI
```
````

### Sink Resolvers

Sink resolver operators may either update the flow graph or the topographic
surface (or both) so that no flow is trapped in closed depressions (i.e., every
flow path is ensured to reach one of the base level nodes).

Some operators like {py:class}`~fastscapelib.PFloodSinkResolver` only update the
topographic surface and do not require any pre-existing flow paths, while other
operators like {py:class}`~fastscapelib.MSTSinkResolver` require single flow
paths that will be updated in-place.

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/flow/sink_resolver.hpp"

namespace fs = fastscapelib;

fs::pflood_sink_resolver::graph_updated      // false
fs::pflood_sink_resolver::elevation_updated  // true
fs::pflood_sink_resolver::in_flowdir         // fs::flow_direction::undefined
fs::pflood_sink_resolver::out_flowdir        // fs::flow_direction::undefined

fs::mst_sink_resolver::graph_updated         // true
fs::mst_sink_resolver::elevation_updated     // true
fs::mst_sink_resolver::in_flowdir            // fs::flow_direction::single
fs::mst_sink_resolver::out_flowdir.          // fs::flow_direction::single
```

```{code-block} Python
:linenos:

import fastscapelib as fs

fs.PFloodSinkResolver.graph_updated       # False
fs.PFloodSinkResolver.elevation_updated   # True
fs.PFloodSinkResolver.in_flowdir          # fs.FlowDirection.UNDEFINED
fs.PFloodSinkResolver.out_flowdir         # fs.FlowDirection.UNDEFINED

fs.MSTSinkResolver.graph_updated       # True
fs.MSTSinkResolver.elevation_updated   # True
fs.MSTSinkResolver.in_flowdir          # fs.FlowDirection.SINGLE
fs.MSTSinkResolver.out_flowdir         # fs.FlowDirection.SINGLE
```
````

### Flow Snapshots

## Flow Routing Strategy (Examples)


## Flow Accumulation

## Advanced Usage

 TODO (flow graph implementation).
