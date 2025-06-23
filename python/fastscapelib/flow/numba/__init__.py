from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Iterable

from _fastscapelib_py.flow import FlowGraphTraversalDir  # type: ignore[import]

if TYPE_CHECKING:
    import numba as nb

    from fastscapelib.flow import FlowGraph, FlowGraphTraversalDir
    from fastscapelib.flow.numba.flow_kernel import (
        NumbaFlowKernel,
        NumbaFlowKernelData,
        NumbaJittedClass,
    )


def create_flow_kernel(
    flow_graph: FlowGraph,
    kernel_func: Callable[[NumbaJittedClass | Any], int | None],
    spec: dict[str, nb.types.Type | tuple[nb.types.Type, Any]],
    *,
    apply_dir: FlowGraphTraversalDir = FlowGraphTraversalDir.DEPTH_UPSTREAM,
    outputs: Iterable[str] = (),
    max_receivers: int = -1,
    n_threads: int = 1,
    print_stats: bool = False,
) -> tuple[NumbaFlowKernel, NumbaFlowKernelData]:
    """Creates a numba flow kernel.

    Parameters
    ----------
    flow_graph : :py:class:`~fastscapelib.FlowGraph`
        A flow graph object.
    kernel_func : callable
        A Python function to apply to each node of the graph. It must take one
        input argument (numba jit-compiled class) that holds or references input
        and output kernel data of one graph node. It should return an integer.
    spec : dict
        Dictionary where keys are kernel input and output variable names and
        values are either variable types (i.e. numba scalar or array types) or
        variable (type, value) tuples (only for scalar variables).
    apply_dir : :py:class:`~fastscapelib.FlowGraphTraversalDir`, optional
        The direction and order in which the flow kernel will be applied along
        the graph (default: downstream to upstream, depth-first search).
    outputs : iterable, optional
        The names of the kernel output variables. All names given here must be
        also present in ``spec``.
    max_receivers : int
        Maximum number of flow receiver nodes per graph node. Setting this number
        to 1 may speed up the application of the kernel function along the graph.
        Setting this number to -1 (default) will use the maximum number defined
        from the flow graph object and is generally safer.
    n_threads : int
        Number of threads to use for applying the kernel function in parallel
        along the flow graph (default: 1, the kernel will be applied sequentially).
        A value > 1 may be useful for kernel functions that are computationally
        expensive. For trivial kernel functions it is generally more efficient to
        apply the kernel sequentially (i.e., ``n_threads = 1``).
    print_stats : bool
        If True, prints a small report on kernel creation performance
        (default: False).

    Returns
    -------
    flow_kernel : :py:class:`~fastscapelib.flow.numba.flow_kernel.NumbaFlowKernel`
        An object used to apply the flow kernel.
    flow_kernel_data : :py:class:`~fastscapelib.flow.numba.flow_kernel.NumbaFlowKernelData`
        An object used to manage flow kernel data.

    """
    from fastscapelib.flow.numba.flow_kernel import NumbaFlowKernelFactory

    factory = NumbaFlowKernelFactory(
        flow_graph,
        kernel_func,
        spec,
        apply_dir,
        outputs=outputs,
        max_receivers=max_receivers,
        n_threads=n_threads,
        print_stats=print_stats,
    )
    return factory.kernel, factory.data
