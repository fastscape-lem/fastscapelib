import pytest

from fastscapelib.flow import (
    DummyFlowRouter,
    FlowGraph,
    NoSinkResolver,
    SingleFlowRouter,
)
from fastscapelib.grid import NodeStatus, ProfileGrid


class TestNoSinkResolver:
    def test___init__(self):
        NoSinkResolver()

        with pytest.raises(TypeError):
            NoSinkResolver(1.0)

    def test_receivers(self):
        profile_grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        FlowGraph(profile_grid, SingleFlowRouter(), NoSinkResolver())
