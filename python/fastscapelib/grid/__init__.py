from typing import TYPE_CHECKING, Any

from _fastscapelib_py.grid import (
    Neighbor,
    Node,
    NodeStatus,
    ProfileBoundaryStatus,
    ProfileGrid,
    RasterBoundaryStatus,
    RasterGrid,
    RasterNeighbor,
    RasterNode,
    TriMesh,
)

try:
    from _fastscapelib_py.grid import HealpixGrid
except ImportError:
    if TYPE_CHECKING:
        raise
    else:
        # Python module built with no healpix grid support (e.g., Windows)
        HealpixGrid: Any = None

__all__ = [
    "HealpixGrid",
    "Neighbor",
    "Node",
    "NodeStatus",
    "ProfileBoundaryStatus",
    "ProfileGrid",
    "RasterBoundaryStatus",
    "RasterGrid",
    "RasterNeighbor",
    "RasterNode",
    "TriMesh",
]
