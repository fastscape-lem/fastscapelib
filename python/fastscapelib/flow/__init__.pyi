import ast
from dataclasses import dataclass
from typing import Any, Callable, ClassVar, Iterable, List, Union, overload

import numba as nb
import numpy as np
import numpy.typing as npt
from _typeshed import Incomplete

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
    def __repr__(self) -> str: ...

class KernelApplicationOrder:
    __members__: ClassVar[dict] = ...  # read-only
    ANY: ClassVar[KernelApplicationOrder] = ...
    DEPTH_DOWNSTREAM: ClassVar[KernelApplicationOrder] = ...
    DEPTH_UPSTREAM: ClassVar[KernelApplicationOrder] = ...
    BREATH_DOWNSTREAM: ClassVar[KernelApplicationOrder] = ...
    BREATH_UPSTREAM: ClassVar[KernelApplicationOrder] = ...
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

class FlowKernelNodeData:
    pass

class FlowKernelData:
    pass

class ConstantAssignmentVisitor(ast.NodeVisitor):
    obj_name: Incomplete
    member_name: Incomplete
    assigned: bool
    aliases: Incomplete
    def __init__(self, obj_name, member_name) -> None: ...
    def visit_Assign(self, node: ast.Assign) -> None: ...
    def visit_AugAssign(self, node: ast.AugAssign) -> None: ...

def timer(msg: str, do_print: bool) -> None: ...

class NumbaKernelData:
    def __init__(
        self,
        grid_size: int,
        spec_keys: list[str],
        grid_data_ty: dict[str, nb.core.types.Type],
        data: FlowKernelData,
    ) -> None: ...
    def __getattr__(self, name): ...
    def __getitem__(self, name): ...
    def __setattr__(self, name, value) -> None: ...
    def __setitem__(self, name, value) -> None: ...
    @property
    def jitclass(self): ...
    @property
    def jitclass_ptr(self): ...
    def bind(self, **kwargs) -> None: ...
    @property
    def bound(self): ...
    def check_bindings(self) -> None: ...

@dataclass
class NumbaKernel:
    kernel: Kernel
    node_data_create: None
    node_data_init: None
    node_data_getter: None
    node_data_setter: None
    func: None
    def __init__(
        self,
        kernel,
        node_data_create,
        node_data_init,
        node_data_getter,
        node_data_setter,
        func,
    ) -> None: ...

def create_flow_kernel(*args, **kwargs) -> tuple[NumbaKernel, NumbaKernelData]: ...

class NumbaFlowKernelFactory:
    def __init__(
        self,
        flow_graph: FlowGraph,
        kernel_func: Callable[[FlowKernelNodeData], int],
        spec: dict[str, nb.core.types.Type | tuple[nb.core.types.Type, Any]],
        application_order: KernelApplicationOrder,
        outputs: Iterable[str] = (),
        max_receivers: int = -1,
        n_threads: int = 1,
        print_generated_code: bool = False,
        print_stats: bool = False,
    ) -> None: ...
    @property
    def kernel(self): ...
    @property
    def data(self): ...
    node_data_setter_tmpl: Incomplete

def py_apply_kernel_impl(
    indices: np.ndarray,
    func: Callable[[FlowKernelNodeData], int],
    data: FlowKernelData,
    node_data: FlowKernelNodeData,
    node_data_getter: Callable[[int, FlowKernelData, FlowKernelNodeData], int],
    node_data_setter: Callable[[int, FlowKernelNodeData, FlowKernelData], int],
): ...
def py_apply_kernel(
    flow_graph: FlowGraph, nb_kernel: NumbaKernel, data: NumbaKernelData
): ...
