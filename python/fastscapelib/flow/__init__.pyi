from typing import Any, ClassVar, List, Union, overload

import numpy as np
import numpy.typing as npt

from fastscapelib.grid import ProfileGrid, RasterGrid, TriMesh

Grid = Union[ProfileGrid, RasterGrid, TriMesh]

class FlowGraphImpl:
    @property
    def single_flow(self) -> bool: ...
    @property
    def receivers(self) -> npt.NDArray[np.uint64]: ...
    @property
    def receivers_count(self) -> npt.NDArray[np.uint64]: ...
    @property
    def receivers_distance(self) -> npt.NDArray[np.float64]: ...
    @property
    def receivers_weight(self) -> npt.NDArray[np.float64]: ...
    @property
    def donors(self) -> npt.NDArray[np.uint64]: ...
    @property
    def donors_count(self) -> npt.NDArray[np.uint64]: ...
    @property
    def dfs_indices(self) -> npt.NDArray[np.uint64]: ...
    @property
    def bfs_indices(self) -> npt.NDArray[np.uint64]: ...
    @property
    def bfs_levels(self) -> npt.NDArray[np.uint64]: ...
    @property
    def basins(self) -> npt.NDArray[np.uint64]: ...

class FlowDirection:
    __members__: ClassVar[dict] = ...  # read-only
    UNDEFINED: ClassVar[FlowDirection] = ...
    SINGLE: ClassVar[FlowDirection] = ...
    MULTI: ClassVar[FlowDirection] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class FlowOperator:
    @property
    def name(self) -> str: ...

class SingleFlowRouter(FlowOperator):
    def __init__(self) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    def __repr__(self) -> str: ...

class MultiFlowRouter(FlowOperator):
    def __init__(self, slope_exp: float = 1.0) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    slope_exp: float
    def __repr__(self) -> str: ...

class PFloodSinkResolver(FlowOperator):
    def __init__(self) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    def __repr__(self) -> str: ...

class MSTMethod:
    __members__: ClassVar[dict] = ...  # read-only
    KRUSKAL: ClassVar[MSTMethod] = ...
    BORUVKA: ClassVar[MSTMethod] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class MSTRouteMethod:
    __members__: ClassVar[dict] = ...  # read-only
    BASIC: ClassVar[MSTRouteMethod] = ...
    CARVE: ClassVar[MSTRouteMethod] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class MSTSinkResolver(FlowOperator):
    def __init__(
        self,
        basin_method: MSTMethod = MSTMethod.KRUSKAL,
        route_method: MSTRouteMethod = MSTRouteMethod.CARVE,
    ) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    basin_method: MSTMethod
    route_method: MSTRouteMethod
    def __repr__(self) -> str: ...

class FlowSnapshot(FlowOperator):
    def __init__(
        self, snapshot_name: str, save_graph: bool = True, save_elevation: bool = False
    ) -> None: ...
    graph_updated: ClassVar[bool] = ...
    elevation_updated: ClassVar[bool] = ...
    in_flowdir: ClassVar[FlowDirection] = ...
    out_flowdir: ClassVar[FlowDirection] = ...
    @property
    def snapshot_name(self) -> str: ...
    @property
    def save_graph(self) -> bool: ...
    @property
    def save_elevation(self) -> bool: ...
    def __repr__(self) -> str: ...

class FlowGraph:
    def __init__(self, grid: Grid, operators: List[FlowOperator]) -> None: ...
    @property
    def operators(self) -> List[FlowOperator]: ...
    @property
    def single_flow(self) -> bool: ...
    def impl(self) -> FlowGraphImpl: ...
    @property
    def graph_snapshot_keys(self) -> List[str]: ...
    def graph_snapshot(self, name: str) -> FlowGraph: ...
    @property
    def elevation_snapshot_keys(self) -> List[str]: ...
    def elevation_snapshot(self, name: str) -> npt.NDArray[np.float64]: ...
    def update_routes(
        self, elevation: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]: ...
    @property
    def base_levels(self) -> list[int]: ...
    @base_levels.setter
    def base_levels(self, value: list[int]) -> None: ...
    @property
    def mask(self) -> npt.NDArray[np.bool_]: ...
    @mask.setter
    def mask(self, value: npt.NDArray[np.bool_]) -> None: ...
    @overload
    def accumulate(
        self, acc: npt.NDArray[np.float64], src: npt.NDArray[np.float64]
    ) -> None: ...
    @overload
    def accumulate(self, acc: npt.NDArray[np.float64], src: float) -> None: ...
    @overload
    def accumulate(self, src: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]: ...
    @overload
    def accumulate(self, src: float) -> npt.NDArray[np.float64]: ...
    def basins(self) -> npt.NDArray[np.uint64]: ...
    @property
    def size(self) -> int: ...
    @property
    def grid_shape(self) -> list[int]: ...
    def __repr__(self) -> str: ...

class FlowGraphTraversalDir:
    __members__: ClassVar[dict] = ...  # read-only
    ANY: ClassVar[FlowGraphTraversalDir] = ...
    DEPTH_DOWNSTREAM: ClassVar[FlowGraphTraversalDir] = ...
    DEPTH_UPSTREAM: ClassVar[FlowGraphTraversalDir] = ...
    BREADTH_DOWNSTREAM: ClassVar[FlowGraphTraversalDir] = ...
    BREADTH_UPSTREAM: ClassVar[FlowGraphTraversalDir] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def data(self) -> Any: ...
    def value(self) -> int: ...

class _FlowKernel:
    def __init__(self) -> None: ...
    @property
    def func(self) -> int: ...
    @func.setter
    def func(self, value: int) -> None: ...
    @property
    def node_data_getter(self) -> int: ...
    @node_data_getter.setter
    def node_data_getter(self, value: int) -> None: ...
    @property
    def node_data_setter(self) -> int: ...
    @node_data_setter.setter
    def node_data_setter(self, value: int) -> None: ...
    @property
    def node_data_create(self) -> int: ...
    @node_data_create.setter
    def node_data_create(self, value: int) -> None: ...
    @property
    def node_data_init(self) -> int: ...
    @node_data_init.setter
    def node_data_init(self, value: int) -> None: ...
    @property
    def node_data_free(self) -> int: ...
    @node_data_free.setter
    def node_data_free(self, value: int) -> None: ...
    @property
    def n_threads(self) -> int: ...
    @n_threads.setter
    def n_threads(self, value: int) -> None: ...
    @property
    def min_block_size(self) -> int: ...
    @min_block_size.setter
    def min_block_size(self, value: int) -> None: ...
    @property
    def min_level_size(self) -> int: ...
    @min_level_size.setter
    def min_level_size(self, value: int) -> None: ...
    @property
    def apply_dir(self) -> FlowGraphTraversalDir: ...
    @apply_dir.setter
    def apply_dir(self, value: FlowGraphTraversalDir) -> None: ...

class _JitClass:
    def __init__(self) -> None: ...
    @property
    def meminfo(self) -> int: ...
    @meminfo.setter
    def meminfo(self, value: int): ...
    @property
    def data(self) -> int: ...
    @data.setter
    def data(self, value: int): ...

class _FlowKernelData:
    def __init__(self) -> None: ...
    @property
    def data(self) -> _JitClass: ...
    @data.setter
    def data(self, value: _JitClass): ...
