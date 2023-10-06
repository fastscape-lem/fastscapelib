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
require updating the flow paths repeatedly) and are accessible through a
flexible and user-friendly API.

(guide-flow-graph)=
## Flow Graph

Flow paths are contained in a "flow graph" data structure (more precisely, a
Directed Acyclic Graph or DAG) that connects each grid node to its flow
receiver(s) (neighbors).

The flow graph data structure is distinct although closely related to a grid
({ref}`Figure 1 <fig_grid_vs_graph>`):

- each graph node correspond to a grid node
- graph edges together form a subset of the grid (implicit or explicit) edges

```{figure} _static/fig_grid_vs_graph.svg
:name: fig_grid_vs_graph
:alt: Flow graph on a raster grid
:align: center

Figure 1: Two kinds of flow graphs (left: single flow direction, right: multiple
flow direction) on top of a raster grid with 8-node connectivity. Dots are grid
(and graph) nodes, dots with a different color at the borders are the base level
nodes, thin dashed lines represent the grid connectivity and thick arrows are
the flow paths (graph edges).
```

A {cpp:class}`~fastscapelib::flow_graph` (C++) or
{py:class}`~fastscapelib.FlowGraph` (Python) object can be created from any grid
object (see Section {ref}`guide-grids`), e.g., in the example below from a
raster grid:

````{tab-set-code}
```{code-block} C++
:linenos:
:emphasize-lines: 10

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;

fs::raster_boundary_status boundaries{ fs::node_status::fixed_value };
fs::raster_grid grid({ 100, 100 }, { 200.0, 200.0 }, boundaries);

fs::flow_graph<fs::raster_grid> graph(grid, { fs::single_flow_router() });
```

```{code-block} Python
:linenos:
:emphasize-lines: 5

import fastscapelib as fs

grid = fs.RasterGrid([100, 100], [200.0, 200.0], fs.NodeStatus.FIXED_VALUE)

graph = fs.FlowGraph(grid, [fs.SingleFlowRouter()])
```
````

A sequence of flow operators must be passed as second argument to the flow graph
constructor in order to define the flow routing strategy (see Sections
{ref}`guide-flow-operators` and {ref}`guide-flow-routing-strategies`).

(guide-base-level-nodes)=
### Base Level Nodes

The flow routing algorithms implemented in Fastscapelib require one or more
nodes of the graph to be set as "base level" (or outlet) nodes ({ref}`Figure 1
<fig_grid_vs_graph>`).

These specific nodes may collect input flow but cannot release any output flow
(it could also mean that the flow is leaving the modeled domain or sub-domain,
e.g., surface land water entering the sea).

When creating a new flow graph, all grid nodes with status
{py:attr}`~fastscapelib.NodeStatus.FIXED_VALUE` are set as base
level nodes by default.

It is possible to set alternative base level nodes via the
{cpp:func}`~fastscapelib::flow_graph::set_base_levels` method (C++) or the
{py:attr}`~fastscapelib.FlowGraph.base_levels` property (Python). This can be
done each time prior to {ref}`updating the flow paths
<guide-compute-flow-paths>` in the cases where the base levels are evolving over
time (e.g., inner lake, sea or ocean dynamics, planet craters).

### Mask

A binary mask can be set via the {cpp:func}`~fastscapelib::flow_graph::set_mask`
method (C++) or the {py:attr}`~fastscapelib.FlowGraph.mask` property (Python) to
exclude grid nodes from the computation of the flow graph. By default, all grid
nodes are included in the flow graph (no mask).

Setting a mask is useful if the modeled domain encompass multiple sub-domains
such as land vs. ocean. Those sub-domains may each have their own flow graph
with different {ref}`flow routing strategies <guide-flow-routing-strategies>`,
{ref}`base level nodes <guide-base-level-nodes>` and masks.

Like for the base level nodes, a flow graph's mask can be set or updated each
time prior to {ref}`updating the flow paths <guide-compute-flow-paths>`, e.g.,
according to sea level change over time.

(guide-compute-flow-paths)=
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

fs::raster_boundary_status boundaries{ fs::node_status::fixed_value };
fs::raster_grid grid({ 100, 100 }, { 200.0, 200.0 }, boundaries);

fs::flow_graph<fs::raster_grid> graph(grid, { fs::single_flow_router() });

xt::xarray<double> elevation = xt::random::rand<double>(grid.shape());

// update_routes returns a const reference!
const auto new_elevation = graph.update_routes(elevation);
```

```{code-block} Python
:linenos:
:emphasize-lines: 10

import numpy as np
import fastscapelib as fs

grid = fs.RasterGrid([100, 100], [200.0, 200.0], fs.NodeStatus.FIXED_VALUE)

graph = fs.FlowGraph(grid, [fs.SingleFlowRouter()])

elevation = np.random.uniform(size=grid.shape)

new_elevation = graph.update_routes(elevation)
```
````

{py:meth}`~fastscapelib.FlowGraph.update_routes` returns another elevation
field, which may differ from the input elevation field depending on the applied
flow operators (e.g., some operators fill the closed depressions found in the
input topography).

:::{note}

{py:meth}`~fastscapelib.FlowGraph.update_routes` never updates in-place the
values of the input elevation field.

 :::

(guide-flow-operators)=
## Flow Operators

Flow can be routed over the topographic surface in many different ways ;
choosing one approach over another highly depends on the case studied.
Fastscapelib relies on the concept of "flow operators" that provide a flexible
and convenient solution for implementing simple to advanced {ref}`flow routing
strategies <guide-flow-routing-strategies>`.

A flow operator is a "routing unit" that:

- may read and/or modify in-place the flow graph instance to which it has been attached
- may read and/or update (a copy of) the elevation values passed to the
  {py:meth}`~fastscapelib.FlowGraph.update_routes` method
- may expose its own parameters

There are currently three categories of operators (see the {ref}`C++
<cpp-api-flow-operators>` and {ref}`Python <py-api-flow-operators>` API
reference for a full list of available operators).

### Flow Routers

Flow router operators generally compute new flow paths from scratch and fully
update the flow graph. The updated flow graph has a defined
{py:class}`~fastscapelib.FlowDirection`: either ``SINGLE`` with one unique
receiver per node or ``MULTI`` where the flow is partitioned among multiple
receiver nodes with fixed or variable weights ({ref}`Figure 1
<fig_grid_vs_graph>`).


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
flow path is ensured to reach one of the {ref}`base level nodes
<guide-base-level-nodes>`).

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

{py:class}`~fastscapelib.FlowSnapshot` is a special operator that can be used to
read and save intermediate states of the flow graph and/or topographic elevation
between two other operators.

````{tab-set-code}
```{code-block} C++
:linenos:

#include "fastscapelib/flow/flow_router.hpp"

namespace fs = fastscapelib;

fs::flow_snapshot::graph_updated      // false
fs::flow_snapshot::elevation_updated  // false
```

```{code-block} Python
:linenos:

import fastscapelib as fs

fs.FlowSnapshot.graph_updated       # False
fs.FlowSnapshot.elevation_updated   # False
```
````

(guide-flow-routing-strategies)=
## Flow Routing Strategies

A flow routing strategy is defined by an ordered sequence of flow operators
passed to the {py:class}`~fastscapelib.FlowGraph` constructor. Calling
{py:meth}`~fastscapelib.FlowGraph.update_routes()` will execute all the
operators one after each other.

Here below are a few examples from simple to more advanced strategies (grid
setup is skipped for more clarity, see the {ref}`flow graph full example
above <guide-flow-graph>`).

### Single Direction Flow Routing

The simplest example is to compute the flow paths along the steepest descent. It
only requires the {py:class}`~fastscapelib.SingleFlowRouter` operator.

````{tab-set-code}
```{code-block} C++
fs::flow_graph<fs::raster_grid> single_graph(grid, { fs::single_flow_router() });
```

```{code-block} Python
single_graph = fs.FlowGraph(grid, [fs.SingleFlowRouter()])
```
````

Applied on a raster grid with 8-node connectivity, this is equivalent to the
so-called "D8" algorithm {cite:p}`OCallaghan1984`.

### Flow Routing Across Closed Depressions

The priority flood sink resolver {py:class}`~fastscapelib.PFloodSinkResolver`
{cite:p}`Barnes2014` can be used to fill the closed depressions in the
topography before computing the flow paths.

````{tab-set-code}
```{code-block} C++
#include "fastscapelib/flow/sink_resolver.hpp"

fs::flow_graph<fs::raster_grid> single_graph_nosink(
    grid, { fs::pflood_sink_resolver(), fs::single_flow_router() });
```

```{code-block} Python
single_graph_nosink = fs.FlowGraph(
    grid, [fs.PFloodSinkResolver(), fs.SingleFlowRouter()]
)
```
````

### Flow Routing Across Closed Depression 2

Alternatively, {py:class}`~fastscapelib.MSTSinkResolver`
{cite:p}`Cordonnier2019` can be used after computing the flow paths to re-route
the flow trapped in closed depressions.

````{tab-set-code}
```{code-block} C++
#include "fastscapelib/flow/sink_resolver.hpp"

fs::flow_graph<fs::raster_grid> single_graph_nosink2(
    grid, { fs::single_flow_router(), fs::mst_sink_resolver() });
```

```{code-block} Python
single_graph_nosink2 = fs.FlowGraph(
    grid, [fs.SingleFlowRouter(), fs.MSTSinkResolver()]
)
```
````

Compared to {py:class}`~fastscapelib.PFloodSinkResolver`,
{py:class}`~fastscapelib.MSTSinkResolver` can be much faster to execute when the
number of closed depressions found in the topography is small (often the case in
LEM simulations already after a few time steps if no process is creating new
closed depressions). However, {py:class}`~fastscapelib.MSTSinkResolver` requires
pre-computed single flow paths while
{py:class}`~fastscapelib.PFloodSinkResolver` has no specific requirement.

### Multiple Direction Flow Routing (with Snapshots)

Here is a more advanced example where depression-free multiple direction flow
paths are obtained after applying several operators. Two graph snapshots are also
saved: the single direction flow paths before and after resolving closed
depressions, respectively.

````{tab-set-code}
```{code-block} C++
#include <memory>
#include "fastscapelib/flow/flow_snapshot.hpp"

auto mrouter_ptr = std::make_shared<fs::multi_flow_router>(1.0);

fs::flow_graph<fs::raster_grid> multi_graph(
    grid,
    { fs::single_flow_router(),
      fs::flow_snapshot("single")
      fs::mst_sink_resolver(),
      fs::flow_snapshot("single_nosink"),
      mrouter_ptr });
```

```{code-block} Python
mrouter = fs.MultiFlowRouter(1.0)

multi_graph = fs.FlowGraph(
    grid,
    [
        fs.SingleFlowRouter(),
        fs.FlowSnapshot("single"),
        fs.MSTSinkResolver(),
        fs.FlowSnapshot("single_nosink"),
        mrouter,
    ],
)
```
````

Note that flow operators may be instantiated outside of the flow graph
constructor, like {py:class}`~fastscapelib.MultiFlowRouter` in the example
above. This is useful for changing parameter values between two successive calls
to {py:meth}`~fastscapelib.FlowGraph.update_routes()`, e.g.,

````{tab-set-code}
```{code-block} C++
multi_graph.update_routes(elevation);

mrouter->m_slope_exp = 1.5;

multi_graph.update_routes(elevation);

```

```{code-block} Python
multi_graph.update_routes(elevation)

mrouter.slope_exp = 1.5

multi_graph.update_routes(elevation)
```
````

The flow partitioning method used in {py:class}`~fastscapelib.MultiFlowRouter`
is equivalent to {cite:t}`Quinn1991` method (``slope_exp = 1``) or
{cite:t}`Holmgren1994` method.

Graph snapshots are separate instances of {py:class}`~fastscapelib.FlowGraph`
accessible via the {py:meth}`~fastscapelib.FlowGraph.graph_snapshot` method of
the main graph instance. They are read-only.

````{tab-set-code}
```{code-block} C++
auto snapshot = multi_graph.graph_snapshot("single");

snapshot.update_routes(elevation)  // Error (read-only)!
```

```{code-block} Python
snapshot = multi_graph.graph_snapshot("single")

snapshot.update_routes(elevation)  # Error (read-only)!
```
````

## Flow Accumulation

After {ref}`computing the flow paths <guide-compute-flow-paths>`, a
{py:class}`~fastscapelib.FlowGraph` object is ready to accumulate some locally
produced quantity or flux along the network via the
{py:meth}`~fastscapelib.FlowGraph.accumulate` method. This is handy for
computing a range of internal variables, model outputs and/or diagnostics during
a simulation. Flow accumulation requires a source term, which may be spatially
uniform or variable.

:::{note}

The input source term passed to {py:meth}`~fastscapelib.FlowGraph.accumulate`
must expressed in units per-area. It will be integrated uniformly over the area
surrounded by each grid node given by ``nodes_areas()``.

:::

A few examples:

- *drainage area (or upslope contributing area)*

````{tab-set-code}
```{code-block} C++
auto drainage_area = graph.accumulate(1.0);
```

```{code-block} Python
drainage_area = graph.accumulate(1.0)
```
````

Where the source term is equal to 1 (unit-less) so that ``drainage_area`` has
dimensions {math}`[L^2]` and only results from the accumulation of the node
(cell) areas along the flow paths.

- *water discharge computed from local surface runoff rate*

````{tab-set-code}
```{code-block} C++
xt::xarray<double> runoff_rate = xt::random::rand<double>(grid.shape());

auto discharge = graph.accumulate(runoff_rate);
```

```{code-block} Python
runoff_rate = np.random.uniform(size=grid.shape)

discharge = graph.accumulate(runoff_rate)
```
````

Where ``runoff_rate`` has dimensions {math}`[L/T]` and ``discharge`` has
dimensions {math}`[L^3/T]`.

- *sediment volume computed from vertical (local) erosion*

````{tab-set-code}
```{code-block} C++
double dt = 1e3;
xt::xarray<double> erosion = channel_eroder.erode(elevation, drainage_area, dt);

auto sediment_vol = graph.accumulate(erosion);
```

```{code-block} Python
dt = 1e3
erosion = channel_eroder.erode(elevation, drainage_area, dt)

sediment_vol = graph.accumulate(erosion)
```
````

Where ``erosion`` has dimensions {math}`[L]` and ``sediment_vol`` has
dimensions {math}`[L^3]`.

## Advanced Usage

For advanced use cases it is possible to have read-only access to the graph
internal data via the {py:meth}`~fastscapelib.FlowGraph.impl` method, which
returns the {py:class}`~fastscapelib.FlowGraphImpl` instance attached to the
flow graph.

:::{warning}

Relying on graph internal data for implementing custom logic in 3rd-party code
should be done only when there is no alternative. Always prefer
{py:meth}`~fastscapelib.FlowGraph.accumulate` when possible. Feedback and/or
feature requests are also
[welcome](https://github.com/fastscape-lem/fastscapelib/discussions).

{py:class}`~fastscapelib.FlowGraphImpl` is an implementation detail and
shouldn't be taken as stable API. In C++ this class is under the
`fastscapelib::detail` namespace.

Graph internal data should **never** be updated directly, or at your own risk!
Only flow operators may do this internally. If direct update is possible this is
probably a bug.

:::
