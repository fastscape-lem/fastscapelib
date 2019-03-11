import numpy as np

import fastscapelib


def test_profile_grid():
    es = [fastscapelib.NodeStatus.FIXED_VALUE_BOUNDARY,
          fastscapelib.NodeStatus.FIXED_VALUE_BOUNDARY]
    g = fastscapelib.ProfileGrid(
        10, 2, es, [(5, fastscapelib.NodeStatus.FIXED_VALUE_BOUNDARY)]
    )

    assert g.size == 10
    assert g.spacing == 2.
    np.testing.assert_equal(g.node_status,
                            np.array([1, 0, 0, 0, 0, 1, 0, 0, 0, 1]))
