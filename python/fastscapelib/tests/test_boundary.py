import numpy as np
import pytest

import fastscapelib
from fastscapelib import NodeStatus


@pytest.fixture(scope='session')
def node_status():
    return np.zeros((5, 5), dtype=np.uint8)


def test_create_node_status(node_status):
    actual_node_status = fastscapelib.create_node_status((5, 5))

    np.testing.assert_array_equal(actual_node_status, node_status)
    assert actual_node_status.dtype == node_status.dtype


def test_set_node_status_grid_boundaries(node_status):
    fastscapelib.set_node_status_grid_boundaries(
        node_status,
        NodeStatus.FIXED_VALUE_BOUNDARY,
        NodeStatus.FIXED_VALUE_BOUNDARY,
        NodeStatus.FIXED_VALUE_BOUNDARY,
        NodeStatus.FIXED_VALUE_BOUNDARY
    )

    assert np.all(
        node_status[[0, -1], :] == int(NodeStatus.FIXED_VALUE_BOUNDARY)
    )
    assert np.all(
        node_status[:, [0, -1]] == int(NodeStatus.FIXED_VALUE_BOUNDARY)
    )

    with pytest.raises(ValueError):
        fastscapelib.set_node_status_grid_boundaries(
            node_status,
            NodeStatus.LOOPED_BOUNDARY,
            NodeStatus.CORE_NODE,
            NodeStatus.CORE_NODE,
            NodeStatus.CORE_NODE
        )
