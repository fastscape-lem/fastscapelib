from typing import ClassVar, List, Tuple, Type, overload

import numpy as np
import numpy.typing as npt

class NodeStatus:
    __members__: ClassVar[dict] = ...  # read-only
    CORE: ClassVar[NodeStatus] = ...
    FIXED_VALUE: ClassVar[NodeStatus] = ...
    FIXED_GRADIENT: ClassVar[NodeStatus] = ...
    LOOPED: ClassVar[NodeStatus] = ...
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

class Node:
    def __init__(self, idx: int, status: NodeStatus) -> None: ...
    idx: int
    status: NodeStatus

class Neighbor:
    def __init__(self, idx: int, distance: float, status: NodeStatus) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    idx: int
    distance: float
    status: NodeStatus

class ProfileBoundaryStatus:
    @overload
    def __init__(self, status: NodeStatus) -> None: ...
    @overload
    def __init__(self, left_status: NodeStatus, right_status: NodeStatus) -> None: ...
    @overload
    def __init__(self, status: List[NodeStatus]) -> None: ...
    left: NodeStatus
    right: NodeStatus
    @property
    def is_horizontal_looped(self) -> bool: ...

class ProfileGrid:
    @overload
    def __init__(
        self,
        size: int,
        spacing: float,
        bounds_status: ProfileBoundaryStatus,
        nodes_status: List[Node],
    ) -> None: ...
    @overload
    def __init__(
        self,
        size: int,
        spacing: float,
        bounds_status: List[NodeStatus],
        nodes_status: List[Tuple[int, NodeStatus]],
    ) -> None: ...
    @classmethod
    def from_length(
        cls: Type[ProfileGrid],
        size: int,
        length: float,
        bounds_status: ProfileBoundaryStatus,
        nodes_status: List[Node],
    ) -> ProfileGrid: ...
    is_structured: ClassVar[bool]
    is_uniform: ClassVar[bool]
    n_neighbors_max: ClassVar[int]
    @property
    def size(self) -> int: ...
    @property
    def shape(self) -> List[int]: ...
    @property
    def spacing(self) -> float: ...
    @property
    def length(self) -> float: ...
    @property
    def nodes_status(self) -> npt.NDArray[np.uint8]: ...
    def neighbors_count(self, idx: int) -> int: ...
    def neighbors_indices(self, idx: int) -> npt.NDArray[np.uint64]: ...
    def neighbors_distances(self, idx: int) -> npt.NDArray[np.float64]: ...
    def neighbors(self, idx: int) -> List[Neighbor]: ...

class RasterNode:
    def __init__(self, row: int, col: int, status: NodeStatus) -> None: ...
    row: int
    col: int
    status: NodeStatus

class RasterNeighbor:
    def __init__(
        self, flatten_idx: int, row: int, col: int, distance: float, status: NodeStatus
    ) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    flatten_idx: int
    row: int
    col: int
    distance: float
    status: NodeStatus

class RasterBoundaryStatus:
    @overload
    def __init__(self, status: NodeStatus) -> None: ...
    @overload
    def __init__(self, status: List[NodeStatus]) -> None: ...
    left: NodeStatus
    right: NodeStatus
    bottom: NodeStatus
    top: NodeStatus
    @property
    def is_horizontal_looped(self) -> bool: ...
    @property
    def is_vertical_looped(self) -> bool: ...

class RasterGrid:
    def __init__(
        self,
        shape: List[int],
        spacing: npt.ArrayLike,
        bounds_status: RasterBoundaryStatus,
        nodes_status: List[RasterNode],
    ) -> None: ...
    @classmethod
    def from_length(
        cls: Type[RasterGrid],
        shape: List[int],
        length: npt.ArrayLike,
        bounds_status: RasterBoundaryStatus,
        nodes_status: List[RasterNode],
    ) -> RasterGrid: ...
    is_structured: ClassVar[bool]
    is_uniform: ClassVar[bool]
    n_neighbors_max: ClassVar[int]
    @property
    def size(self) -> int: ...
    @property
    def shape(self) -> List[int]: ...
    @property
    def spacing(self) -> float: ...
    @property
    def length(self) -> float: ...
    @property
    def nodes_status(self) -> npt.NDArray[np.uint8]: ...
    def neighbors_count(self, idx: int) -> int: ...
    @overload
    def neighbors_indices(self, idx: int) -> npt.NDArray[np.uint64]: ...
    @overload
    def neighbors_indices(self, row: int, col: int) -> npt.NDArray[np.uint64]: ...
    def neighbors_distances(self, idx: int) -> npt.NDArray[np.float64]: ...
    @overload
    def neighbors(self, idx: int) -> List[Neighbor]: ...
    @overload
    def neighbors(self, row: int, col: int) -> List[RasterNeighbor]: ...

class UnstructuredMesh:
    def __init__(
        self,
        points: npt.NDArray[np.float64],
        neighbors_indices_ptr: npt.NDArray[np.uint64],
        neighbors_indices: npt.NDArray[np.uint64],
        convex_hull_indices: npt.NDArray[np.uint64],
        areas: npt.NDArray[np.float64],
        nodes_status: List[Node],
    ) -> None: ...
    is_structured: ClassVar[bool]
    is_uniform: ClassVar[bool]
    n_neighbors_max: ClassVar[int]
    @property
    def size(self) -> int: ...
    @property
    def shape(self) -> List[int]: ...
    @property
    def nodes_status(self) -> npt.NDArray[np.uint8]: ...
    def neighbors_count(self, idx: int) -> int: ...
    def neighbors_indices(self, idx: int) -> npt.NDArray[np.uint64]: ...
    def neighbors_distances(self, idx: int) -> npt.NDArray[np.float64]: ...
    def neighbors(self, idx: int) -> List[Neighbor]: ...
