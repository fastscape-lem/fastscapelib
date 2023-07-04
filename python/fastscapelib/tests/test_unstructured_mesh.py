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
    triangles = np.array([[2, 4, 3], [4, 2, 1], [4, 0, 3], [0, 4, 1]])

    # for simplicity, let's set the area of boundary nodes all equal to 1
    # and the area of the inner node equal to 2
    areas = np.array([1.0, 1.0, 1.0, 1.0, 2.0])

    return {
        "points": points,
        "triangles": triangles,
        "areas": areas,
    }


class TestUnstructuredMesh:
    def test_static_properties(self) -> None:
        assert UnstructuredMesh.is_structured is False
        assert UnstructuredMesh.is_uniform is False
        assert UnstructuredMesh.n_neighbors_max == 20

    def test_constructor(self, mesh_args) -> None:
        mesh = UnstructuredMesh(*mesh_args.values(), [])  # type: ignore[call-arg]

        assert mesh.size == 5
        assert mesh.shape == [5]

    def test_nodes_indices(self, mesh_args) -> None:
        mesh = UnstructuredMesh(*mesh_args.values(), [])  # type: ignore[call-arg]
        npt.assert_equal(mesh.nodes_indices(), np.arange(mesh.size))
        npt.assert_equal(mesh.nodes_indices(NodeStatus.FIXED_VALUE), [0, 1, 2, 3])
        npt.assert_equal(mesh.nodes_indices(NodeStatus.CORE), [4])
        assert not len(mesh.nodes_indices(NodeStatus.FIXED_GRADIENT))

    def test_nodes_status_default(self, mesh_args) -> None:
        # all boundary nodes (convex hull) have fixed value status
        mesh = UnstructuredMesh(*mesh_args.values(), [])  # type: ignore[call-arg]

        actual = mesh.nodes_status()
        expected = [1, 1, 1, 1, 0]

        npt.assert_array_equal(actual, expected)

        # should return a copy
        assert not np.shares_memory(mesh.nodes_status(), mesh.nodes_status())

        assert mesh.nodes_status(0) == NodeStatus.FIXED_VALUE
        assert mesh.nodes_status(4) == NodeStatus.CORE

    def test_nodes_status_custom(self, mesh_args) -> None:
        bc = [Node(2, NodeStatus.FIXED_VALUE)]
        mesh = UnstructuredMesh(*mesh_args.values(), bc)  # type: ignore[call-arg]

        actual = mesh.nodes_status()
        expected = np.zeros(mesh.size, dtype=np.uint8)
        expected[2] = 1

        npt.assert_array_equal(actual, expected)

        with pytest.raises(ValueError, match=".*not allowed.*"):
            UnstructuredMesh(*mesh_args.values(), [Node(2, NodeStatus.LOOPED)])  # type: ignore[call-arg]

    def test_nodes_areas(self, mesh_args) -> None:
        mesh = UnstructuredMesh(*mesh_args.values(), [])  # type: ignore[call-arg]

        assert mesh.nodes_areas(0) == mesh_args["areas"][0]
        assert mesh.nodes_areas(4) == mesh_args["areas"][4]
        npt.assert_equal(mesh.nodes_areas(), mesh_args["areas"])

    def test_neighbors_count(self, mesh_args) -> None:
        mesh = UnstructuredMesh(*mesh_args.values(), [])  # type: ignore[call-arg]

        assert mesh.neighbors_count(4) == 4
        assert mesh.neighbors_count(0) == 3

    def test_neighbors_indices(self, mesh_args) -> None:
        mesh = UnstructuredMesh(*mesh_args.values(), [])  # type: ignore[call-arg]

        npt.assert_equal(np.sort(mesh.neighbors_indices(3)), [0, 2, 4])
        npt.assert_equal(np.sort(mesh.neighbors_indices(4)), [0, 1, 2, 3])

    def test_neighbors_distances(self, mesh_args) -> None:
        mesh = UnstructuredMesh(*mesh_args.values(), [])  # type: ignore[call-arg]

        dist_diag = np.sqrt(0.5**2 + 0.5**2)

        sorted_nidx = np.argsort(mesh.neighbors_indices(3)).ravel()
        npt.assert_equal(
            mesh.neighbors_distances(3)[sorted_nidx], [dist_diag, dist_diag, 0.5]
        )

        npt.assert_equal(mesh.neighbors_distances(4), [0.5] * 4)

    def test_neighbors(self, mesh_args) -> None:
        mesh = UnstructuredMesh(*mesh_args.values(), [])  # type: ignore[call-arg]

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
