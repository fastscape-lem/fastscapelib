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
    n_threads: int = 1,
    get_data_at_receivers: bool = True,
    set_data_at_receivers: bool = True,
    get_data_at_donors: bool = True,
    set_data_at_donors: bool = True,
    auto_resize: bool = False,
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
        and output kernel data of one graph node. It may return an integer.
    spec : dict
        Dictionary where keys are kernel input and output variable names and
        values are either variable types (i.e. numba scalar or array types) or
        variable (type, value) tuples (only for scalar variables).
    apply_dir : :py:class:`~fastscapelib.FlowGraphTraversalDir`, optional
        The direction and order in which the flow kernel will be applied along
        the graph (default: downstream to upstream, depth-first search).
    outputs : iterable
        The names of the kernel output variables. All names given here must be
        also present in ``spec``. There must be at least one output.
    n_threads : int, optional
        Number of threads to use for applying the kernel function in parallel
        along the flow graph (default: 1, the kernel will be applied sequentially).
        A value > 1 may be useful for kernel functions that are computationally
        expensive. For trivial kernel functions it is generally more efficient to
        apply the kernel sequentially (i.e., ``n_threads = 1``).
    get_data_at_receivers : bool, optional
        If True (default), kernel input and output (array) variables will
        have their values accessible at the node ``.receivers`` from within
        the kernel function (i.e., values are copied from the kernel data prior to
        calling the kernel function). You can set it to False if this is not needed,
        it may speed-up the application of the kernel.
    set_data_at_receivers : bool, optional
        If True (default), allow setting the values of output variables
        at the node ``.receivers`` from within the kernel function (i.e., output values
        are copied back to the kernel data after calling the kernel function).
        You can set it to False if this is not needed, it may speed-up the application
        of the kernel.
        Note: when the kernel function is evaluated at a base level node, any value set
        via ``.receivers`` within the function won't be copied in the output variable.
    get_data_at_donors : bool, optional
        Same than ``get_data_at_receivers`` but for the node ``.donors``.
    set_data_at_donors : bool, optional
        Same than ``set_data_at_receivers`` but for the node ``.donors``.
    auto_resize: bool, optional
        If True, dynamically resize the temporary containers used to store data at flow
        receivers and donors. Otherwise (default), set fixed sizes according to the
        flow graph properties. Default settings are generally recommended.
    print_stats : bool
        If True, prints a small report on kernel creation performance
        (default: False).

    Returns
    -------
    flow_kernel : :py:class:`~fastscapelib.flow.numba.flow_kernel.NumbaFlowKernel`
        Flow kernel proxy object.
    flow_kernel_data : :py:class:`~fastscapelib.flow.numba.flow_kernel.NumbaFlowKernelData`
        Flow kernel data management.

    """
    from fastscapelib.flow.numba.flow_kernel import NumbaFlowKernelFactory

    factory = NumbaFlowKernelFactory(
        flow_graph,
        kernel_func,
        spec,
        apply_dir=apply_dir,
        outputs=outputs,
        n_threads=n_threads,
        get_data_at_receivers=get_data_at_receivers,
        set_data_at_receivers=set_data_at_receivers,
        get_data_at_donors=get_data_at_donors,
        set_data_at_donors=set_data_at_donors,
        auto_resize=auto_resize,
        print_stats=print_stats,
    )
    return factory.kernel, factory.data
