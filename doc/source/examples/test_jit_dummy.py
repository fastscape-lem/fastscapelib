import fastscapelib as fs
import matplotlib.pyplot as plt
import numpy as np

import time
import numba as nb
from pyfs import NumbaFlowKernel, py_apply_kernel
from testpy import Kernel, GraphRunner

grid = fs.RasterGrid.from_length([100, 1000], [5e4, 7.5e4], fs.NodeStatus.FIXED_VALUE)
flow_graph = fs.FlowGraph(grid, [fs.SingleFlowRouter(), fs.MSTSinkResolver()])
rng = np.random.Generator(np.random.PCG64(1234))
foo = rng.uniform(0, 5, size=grid.shape).ravel()
flow_graph.update_routes(foo)

# print(flow_graph.impl().bfs_indices)
# print(flow_graph.impl().bfs_levels)

# print(flow_graph.impl().bfs_levels[1:] - flow_graph.impl().bfs_levels[:-1])

# print(flow_graph.impl().bfs_levels.shape)
# print(np.all(flow_graph.impl().receivers_count == 1))


dt = 2e4


def kernel_func(node, dt):
    node.foo += 1.0
    node.baz += 1

    for i in range(100):
        node.bar *= node.foo * np.power(1.1, node.foo * i) * np.exp(node.foo * i)


foo = np.ones(flow_graph.size)
baz = np.zeros(flow_graph.size, dtype=np.int64)

kernel = NumbaFlowKernel(
    flow_graph,
    kernel_func,
    grid_data={
        "foo": foo,
        "bar": np.zeros(flow_graph.size),
        "baz": baz,
    },
    constants={},
    outputs=["foo", "baz"],
    max_receivers=1,
    n_threads=1,
    # print_generated_code=True,
)

N_TIMES = 1

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
