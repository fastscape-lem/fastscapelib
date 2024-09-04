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

__all__ = [
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

try:
    from _fastscapelib_py import HealpixGrid

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
except ImportError:
    # Python module built with no healpix grid support (e.g., Windows)
    pass
