import fastscapelib as fs
import matplotlib.pyplot as plt
import numpy as np

import numba as nb
from pyfs import NumbaFlowKernel, py_apply_kernel
from testpy import GraphRunner

grid = fs.RasterGrid.from_length([1000, 1000], [5e4, 7.5e4], fs.NodeStatus.FIXED_VALUE)
flow_graph = fs.FlowGraph(grid, [fs.SingleFlowRouter(), fs.MSTSinkResolver()])

rng = np.random.Generator(np.random.PCG64(1234))

init_elevation = rng.uniform(0, 5, size=grid.shape)
drainage_area = np.empty_like(init_elevation)

uplift_rate = np.full_like(init_elevation, 1e-3)
uplift_rate[[0, -1], :] = 0.0
uplift_rate[:, [0, -1]] = 0.0

nsteps = 50
dt = 2e4
flow_graph.update_routes(init_elevation)

from numba import float64, int32, int64, njit, typeof, uint64
from numba.experimental import jitclass


def kernel_func(node, dt):

    r_count = node.receivers.count
    if r_count == 1 and node.receivers.distance[0] == 0.0:
        return

    elevation_flooded = np.finfo(np.double).max

    for r in range(r_count):
        irec_elevation_next = node.receivers.elevation[r] - node.receivers.erosion[r]

        if irec_elevation_next < elevation_flooded:
            elevation_flooded = irec_elevation_next

    if node.elevation <= elevation_flooded:
        return

    eq_num = node.elevation
    eq_den = 1.0

    for r in range(r_count):
        irec_elevation = node.receivers.elevation[r]
        irec_elevation_next = irec_elevation - node.receivers.erosion[r]

        if irec_elevation > node.elevation:
            continue

        irec_weight = node.receivers.weight[r]
        irec_distance = node.receivers.distance[r]

        factor = (
            node.k_coef * dt * np.power(node.drainage_area * irec_weight, node.area_exp)
        )
        factor /= irec_distance
        eq_num += factor * irec_elevation_next
        eq_den += factor

    elevation_updated = eq_num / eq_den

    if elevation_updated < elevation_flooded:
        elevation_updated = elevation_flooded + np.finfo(np.double).tiny

    node.erosion = node.elevation - elevation_updated


kernel = NumbaFlowKernel(
    flow_graph,
    kernel_func,
    grid_data={
        "drainage_area": np.ones(flow_graph.size),
        "erosion": np.zeros(flow_graph.size),
        "elevation": init_elevation.ravel().copy(),
    },
    constants={"k_coef": 2e-4, "area_exp": 0.4, "slope_exp": 1.0},
    outputs=["erosion"],
    max_receivers=1,
    n_threads=1,
    # print_generated_code=True,
    # print_stats=True,
)

import time

N_TIMES = 100
indices = flow_graph.impl().bfs_indices
levels = flow_graph.impl().bfs_levels
print("levels:", levels.size)
blocks = levels[1:] - levels[:-1]
print(
    "min/max/median/mean levels size:",
    blocks.min(),
    blocks.max(),
    round(np.median(blocks)),
    round(blocks.mean()),
)


def bench_fn(fn, kernel):

    n_threads = kernel.kernel.n_threads
    if fn != py_apply_kernel:
        kernel = kernel.kernel
        msg = "C++ apply kernel"
    else:
        msg = "Python apply_kernel"

    t1 = time.time()
    for _ in range(N_TIMES):
        fn(indices, levels, kernel, dt)
    elapsed_time = time.time() - t1
    print(
        f"{msg} ({n_threads} threads) {round(elapsed_time, 2)}s,",
        f"{round(elapsed_time * 1e9 / N_TIMES / grid.size)}ns/op,",
        f"{round(elapsed_time * 1e9 / N_TIMES / grid.size * n_threads)}ns/op/thread",
    )


bench_fn(py_apply_kernel, kernel)

runner = GraphRunner(1)
bench_fn(runner.apply_kernel, kernel)
del runner

for n_threads in [2, 4, 8, 16, 22]:
    runner = GraphRunner(n_threads)
    kernel.kernel.n_threads = n_threads
    bench_fn(
        runner.apply_kernel,
        kernel,
    )
    del runner
