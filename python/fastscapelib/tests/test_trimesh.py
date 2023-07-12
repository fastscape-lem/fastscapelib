import numpy as np
import numpy.testing as npt
import pytest

from fastscapelib.grid import Neighbor, NodeStatus, TriMesh


@pytest.fixture
def mesh_args():
    # simple mesh with four boundary nodes and one inner node:
    #      o
    #   o  o  o
    #      o
    #
    # Returns arguments passed to the TriMesh constructor

    points = np.array([[0.0, 0.5], [0.5, 0.0], [0.0, -0.5], [-0.5, 0.0], [0.0, 0.0]])
    triangles = np.array([[2, 4, 3], [4, 2, 1], [4, 0, 3], [0, 4, 1]])

    return points, triangles


class TestTriMesh:
    def test_static_properties(self) -> None:
        assert TriMesh.is_structured is False
        assert TriMesh.is_uniform is False
        assert TriMesh.n_neighbors_max == 20

    def test_constructor(self, mesh_args) -> None:
        points, triangles = mesh_args
        mesh = TriMesh(points, triangles)

        assert mesh.size == 5
        assert mesh.shape == [5]

        with pytest.raises(ValueError, match="invalid shape for points array"):
            TriMesh(np.array([[0.0, 0.5, 1.0]]), triangles)

        with pytest.raises(ValueError, match="invalid shape for triangles array"):
            TriMesh(points, np.array([[0, 1]]))

        with pytest.raises(ValueError, match="invalid shape for nodes_status array"):
            TriMesh(points, triangles, [0])

    def test_nodes_indices(self, mesh_args) -> None:
        mesh = TriMesh(*mesh_args)
        npt.assert_equal(mesh.nodes_indices(), np.arange(mesh.size))
        npt.assert_equal(mesh.nodes_indices(NodeStatus.FIXED_VALUE), [0, 1, 2, 3])
        npt.assert_equal(mesh.nodes_indices(NodeStatus.CORE), [4])
        assert not len(mesh.nodes_indices(NodeStatus.FIXED_GRADIENT))

    def test_nodes_status_default(self, mesh_args) -> None:
        points, triangles = mesh_args

        # all boundary nodes have fixed value status
        mesh = TriMesh(points, triangles)

        actual = mesh.nodes_status()
        expected = [1, 1, 1, 1, 0]

        npt.assert_array_equal(actual, expected)

        # should return a copy
        assert not np.shares_memory(mesh.nodes_status(), mesh.nodes_status())

        assert mesh.nodes_status(0) == NodeStatus.FIXED_VALUE
        assert mesh.nodes_status(4) == NodeStatus.CORE

        # initialize with an array of node status
        mesh2 = TriMesh(points, triangles, [1, 1, 0, 1, 1])
        npt.assert_array_equal(mesh2.nodes_status(), [1, 1, 0, 1, 1])

    def test_nodes_status_custom(self, mesh_args) -> None:
        points, triangles = mesh_args
        mesh = TriMesh(points, triangles, {2: NodeStatus.FIXED_VALUE})

        actual = mesh.nodes_status()
        expected = np.zeros(mesh.size, dtype=np.uint8)
        expected[2] = 1

        npt.assert_array_equal(actual, expected)

        with pytest.raises(ValueError, match=".*not allowed.*"):
            TriMesh(points, triangles, {2: NodeStatus.LOOPED})

    def test_nodes_areas(self) -> None:
        # simple case of 1x1 domain with 9 evenly spaced nodes
        points = np.array(
            [
                [0, 0],
                [0, 0.5],
                [0, 1],
                [0.5, 1],
                [1, 1],
                [1, 0.5],
                [1, 0],
                [0.5, 0],
                [0.5, 0.5],
            ]
        )
        triangles = np.array(
            [
                [0, 1, 8],
                [1, 2, 8],
                [2, 3, 8],
                [3, 4, 8],
                [4, 5, 8],
                [5, 6, 8],
                [6, 7, 8],
                [7, 0, 8],
            ]
        )
        mesh = TriMesh(points, triangles)

        assert mesh.nodes_areas(0) == 0.0625
        assert mesh.nodes_areas(8) == 0.25
        npt.assert_equal(
            mesh.nodes_areas(),
            [0.0625, 0.125, 0.0625, 0.125, 0.0625, 0.125, 0.0625, 0.125, 0.25],
        )

    def test_neighbors_count(self, mesh_args) -> None:
        mesh = TriMesh(*mesh_args)

        assert mesh.neighbors_count(4) == 4
        assert mesh.neighbors_count(0) == 3

    def test_neighbors_indices(self, mesh_args) -> None:
        mesh = TriMesh(*mesh_args)

        npt.assert_equal(np.sort(mesh.neighbors_indices(3)), [0, 2, 4])
        npt.assert_equal(np.sort(mesh.neighbors_indices(4)), [0, 1, 2, 3])

    def test_neighbors_distances(self, mesh_args) -> None:
        mesh = TriMesh(*mesh_args)

        dist_diag = np.sqrt(0.5**2 + 0.5**2)

        sorted_nidx = np.argsort(mesh.neighbors_indices(3)).ravel()
        npt.assert_equal(
            mesh.neighbors_distances(3)[sorted_nidx], [dist_diag, dist_diag, 0.5]
        )

        npt.assert_equal(mesh.neighbors_distances(4), [0.5] * 4)

    def test_neighbors(self, mesh_args) -> None:
        mesh = TriMesh(*mesh_args)

        dist_diag = np.sqrt(0.5**2 + 0.5**2)

        sorted_nidx = np.argsort(mesh.neighbors_indices(3)).ravel()
        temp = mesh.neighbors(3)
        actual = [temp[i] for i in sorted_nidx]

        expected = [
            Neighbor(0, dist_diag, NodeStatus.FIXED_VALUE),
            Neighbor(2, dist_diag, NodeStatus.FIXED_VALUE),
            Neighbor(4, 0.5, NodeStatus.CORE),
        ]

        assert actual == expected

        with pytest.raises(IndexError, match="grid index out of range"):
            mesh.neighbors(111)
