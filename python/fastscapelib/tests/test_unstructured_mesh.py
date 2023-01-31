import numpy as np
import numpy.testing as npt
import pytest
from fastscapelib.grid import Neighbor, Node, NodeStatus, UnstructuredMesh


@pytest.fixture
def mesh_args():
    # simple mesh with four boundary nodes and one inner node:
    #      o
    #   o  o  o
    #      o
    #
    # Returns arguments passed to the UnstructuredMesh constructor

    points = np.array([[0.0, 0.5], [0.5, 0.0], [0.0, -0.5], [-0.5, 0.0], [0.0, 0.0]])

    # the arrays below have been obtained using scipy's Delaunay class (not imported here
    # to avoid a heavy dependency)
    #
    # from scipy.spatial import Delaunay
    # tri = Delaunay(points)
    # indptr, indices = tri.vertex_neighbor_vertices
    # convex_hull_vertices = np.unique(tri.convex_hull.ravel())
    indptr = np.array([0, 3, 6, 9, 12, 16], dtype=np.int32)
    indices = np.array([4, 3, 1, 4, 2, 0, 4, 3, 1, 2, 4, 0, 2, 3, 1, 0], dtype=np.int32)
    convex_hull_indices = np.array([0, 1, 2, 3], dtype=np.int32)

    # for simplicity, let's set the area of boundary nodes all equal to 1
    # and the area of the inner node equal to 2
    areas = np.array([1.0, 1.0, 1.0, 1.0, 2.0])

    return {
        "points": points,
        "indptr": indptr,
        "indices": indices,
        "convex_hull_indices": convex_hull_indices,
        "areas": areas,
    }


class TestUnstructuredMesh:
    def test_static_properties(self):
        assert UnstructuredMesh.is_structured is False
        assert UnstructuredMesh.is_uniform is False
        assert UnstructuredMesh.n_neighbors_max == 30

    def test_constructor(self, mesh_args):
        mesh = UnstructuredMesh(*mesh_args.values(), [])

        assert mesh.size == 5
        assert mesh.shape == [5]

    def test_status_at_nodes_default(self, mesh_args):
        # all boundary nodes (convex hull) have fixed value status
        mesh = UnstructuredMesh(*mesh_args.values(), [])

        actual = mesh.status_at_nodes
        expected = np.zeros(mesh.size, dtype=np.uint8)
        expected[mesh_args["convex_hull_indices"]] = 1

        npt.assert_array_equal(actual, expected)

    def test_status_at_nodes_custom(self, mesh_args):
        bc = [Node(2, NodeStatus.FIXED_VALUE_BOUNDARY)]
        mesh = UnstructuredMesh(*mesh_args.values(), bc)

        actual = mesh.status_at_nodes
        expected = np.zeros(mesh.size, dtype=np.uint8)
        expected[2] = 1

        npt.assert_array_equal(actual, expected)

        with pytest.raises(ValueError, match=".*not allowed.*"):
            UnstructuredMesh(*mesh_args.values(), [Node(2, NodeStatus.LOOPED_BOUNDARY)])

    def test_neighbors_count(self, mesh_args):
        mesh = UnstructuredMesh(*mesh_args.values(), [])

        assert mesh.neighbors_count(4) == 4
        assert mesh.neighbors_count(0) == 3

    def test_neighbors_indices(self, mesh_args):
        mesh = UnstructuredMesh(*mesh_args.values(), [])

        npt.assert_equal(mesh.neighbors_indices(3), [2, 4, 0])
        npt.assert_equal(mesh.neighbors_indices(4), [2, 3, 1, 0])

    def test_neighbors_distances(self, mesh_args):
        mesh = UnstructuredMesh(*mesh_args.values(), [])

        dist_diag = np.sqrt(0.5**2 + 0.5**2)

        npt.assert_equal(mesh.neighbors_distances(3), [dist_diag, 0.5, dist_diag])
        npt.assert_equal(mesh.neighbors_distances(4), [0.5] * 4)

    def test_neighbors(self, mesh_args):
        mesh = UnstructuredMesh(*mesh_args.values(), [])

        dist_diag = np.sqrt(0.5**2 + 0.5**2)

        assert mesh.neighbors(3) == [
            Neighbor(2, dist_diag, NodeStatus.FIXED_VALUE_BOUNDARY),
            Neighbor(4, 0.5, NodeStatus.CORE),
            Neighbor(0, dist_diag, NodeStatus.FIXED_VALUE_BOUNDARY),
        ]

        with pytest.raises(IndexError, match="grid index out of range"):
            mesh.neighbors(111)
