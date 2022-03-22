import numpy as np
import numpy.testing as npt
import pytest

from fastscapelib.flow import (
    DummyFlowRouter,
    FlowGraph,
    FlowRouterMethods,
    MultipleFlowRouter,
    NoSinkResolver,
    SingleFlowRouter,
    SinkResolverMethods,
)
from fastscapelib.grid import (
    Node,
    NodeStatus,
    ProfileGrid,
    RasterBoundaryStatus,
    RasterGrid,
    RasterNode,
)


class TestSinkResolverMethods:
    def test___init__(self):
        SinkResolverMethods.NONE
        SinkResolverMethods.FILL_PFLOOD
        SinkResolverMethods.FILL_MST_KRUSKAL
        SinkResolverMethods.FILL_MST_BORUVKA
        SinkResolverMethods.FILL_AUTO
        SinkResolverMethods.CARVE_MST_KRUSKAL
        SinkResolverMethods.CARVE_MST_BORUVKA


class TestNoSinkResolver:
    def test___init__(self):
        NoSinkResolver()

        with pytest.raises(TypeError):
            NoSinkResolver(1.0)

    def test_receivers(self):
        profile_grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        flow_graph = FlowGraph(profile_grid, SingleFlowRouter(), NoSinkResolver())
