from textwrap import dedent

import numpy as np
import numpy.testing as npt
import pytest
from fastscapelib.flow import (
    FlowGraph,
    FlowSnapshot,
    MultiFlowRouter,
    PFloodSinkResolver,
    SingleFlowRouter,
)
from fastscapelib.grid import NodeStatus, ProfileGrid, RasterBoundaryStatus, RasterGrid


class TestFlowGraph:
    def test___init__(self) -> None:
        profile_grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        raster_grid = RasterGrid(
            [5, 10],
            [2.2, 2.4],
            RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY),
            [],
        )

        FlowGraph(profile_grid, [SingleFlowRouter()])
        FlowGraph(raster_grid, [PFloodSinkResolver(), SingleFlowRouter()])

        with pytest.raises(
            ValueError, match="must have at least one operator that updates"
        ):
            FlowGraph(raster_grid, [])

        with pytest.raises(TypeError, match="invalid flow operator"):
            FlowGraph(raster_grid, ["not a flow operator"])  # type: ignore[list-item]

    def test_operators(self) -> None:
        grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])

        resolver = PFloodSinkResolver()
        router = SingleFlowRouter()

        flow_graph = FlowGraph(grid, [resolver, router])

        assert flow_graph.operators == [resolver, router]

    def test_single_flow(self) -> None:
        grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])

        flow_graph = FlowGraph(grid, [SingleFlowRouter()])
        assert flow_graph.single_flow is True

        flow_graph = FlowGraph(grid, [MultiFlowRouter(1.0)])
        assert flow_graph.single_flow is False

    def test_repr(self) -> None:
        grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        flow_graph = FlowGraph(grid, [SingleFlowRouter(), FlowSnapshot("test")])

        actual = repr(flow_graph)
        expected = """
        <FlowGraph (8 nodes)>
        Operators:
            SingleFlowRouter
            FlowSnapshot 'test' (graph=True, elevation=False)
        """
        assert actual == dedent(expected[1:])

        actual = repr(flow_graph.graph_snapshot("test"))
        expected = """
        <FlowGraph (8 nodes)>
        Operators:
            *empty*
        """
        assert actual == dedent(expected[1:])

    def test_update_routes(self) -> None:
        grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        flow_graph = FlowGraph(grid, [SingleFlowRouter()])

        # pit at 3rd node
        elevation = np.array([0.0, 0.2, 0.1, 0.2, 0.4, 0.6, 0.3, 0.0])

        new_elevation = flow_graph.update_routes(elevation)

        npt.assert_array_equal(elevation, new_elevation)

        npt.assert_equal(
            flow_graph.impl().receivers[:, 0], np.array([0, 0, 2, 2, 3, 6, 7, 7])
        )
        npt.assert_equal(flow_graph.impl().receivers_count, np.ones(elevation.size))
        npt.assert_equal(
            flow_graph.impl().receivers_weight[:, 0], np.ones(elevation.size)
        )

        m = np.iinfo(np.uint64).max
        npt.assert_equal(
            flow_graph.impl().donors,
            np.array(
                [
                    [1, m, m],
                    [m, m, m],
                    [2, 3, m],
                    [4, m, m],
                    [m, m, m],
                    [m, m, m],
                    [5, m, m],
                    [6, m, m],
                ]
            ),
        )
        npt.assert_equal(
            flow_graph.impl().donors_count, np.array([1, 0, 2, 1, 0, 0, 1, 1])
        )

    def test_accumulate_basins(self) -> None:
        # --- test profile grid
        grid = ProfileGrid(8, 2.0, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        flow_graph = FlowGraph(grid, [SingleFlowRouter()])

        # pit at 3rd node
        elevation = np.array([0.0, 0.2, 0.1, 0.2, 0.4, 0.6, 0.3, 0.0])

        new_elevation = flow_graph.update_routes(elevation)

        npt.assert_array_equal(elevation, new_elevation)

        acc = np.empty_like(elevation)
        src = np.ones_like(elevation)
        expected = np.array([4.0, 2.0, 6.0, 4.0, 2.0, 2.0, 4.0, 6.0])

        flow_graph.accumulate(acc, 1.0)
        npt.assert_array_equal(acc, expected)

        flow_graph.accumulate(acc, src)
        npt.assert_array_equal(acc, expected)

        npt.assert_almost_equal(flow_graph.accumulate(1.0), expected)
        npt.assert_almost_equal(flow_graph.accumulate(src), expected)

        # bottom border base-level
        bottom_base_level = [
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.FIXED_VALUE_BOUNDARY,
        ]

        # --- test raster grid
        rgrid = RasterGrid(
            [4, 4],
            [1.0, 1.0],
            RasterBoundaryStatus(bottom_base_level),
            [],
        )
        flow_graph = FlowGraph(rgrid, [SingleFlowRouter()])

        # planar surface tilted along the y-axis + small carved channel
        elevation = np.array(
            [
                [0.6, 0.6, 0.6, 0.6],
                [0.4, 0.4, 0.4, 0.4],
                [0.2, 0.2, 0.2, 0.2],
                [0.1, 0.0, 0.1, 0.1],
            ]
        )

        new_elevation = flow_graph.update_routes(elevation)
        npt.assert_array_equal(elevation, new_elevation)

        acc = np.empty_like(elevation)
        src = np.ones_like(elevation)
        expected = np.array(
            [
                [1.0, 1.0, 1.0, 1.0],
                [2.0, 2.0, 2.0, 2.0],
                [3.0, 3.0, 3.0, 3.0],
                [1.0, 10.0, 1.0, 4.0],
            ]
        )

        flow_graph.accumulate(acc, 1.0)
        npt.assert_array_equal(acc, expected)

        flow_graph.accumulate(acc, src)
        npt.assert_array_equal(acc, expected)

        npt.assert_almost_equal(flow_graph.accumulate(1.0), expected)
        npt.assert_almost_equal(flow_graph.accumulate(src), expected)

        expected = np.array(
            [
                [1, 1, 1, 3],
                [1, 1, 1, 3],
                [1, 1, 1, 3],
                [0, 1, 2, 3],
            ]
        )
        npt.assert_array_equal(flow_graph.basins(), expected)
        npt.assert_array_equal(flow_graph.basins().flatten(), flow_graph.impl().basins)
