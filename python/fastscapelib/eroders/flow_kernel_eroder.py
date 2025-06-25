from __future__ import annotations

import abc
import inspect
from typing import TYPE_CHECKING, Any

import numpy as np

from fastscapelib.flow import FlowGraph, FlowGraphTraversalDir
from fastscapelib.flow.numba import create_flow_kernel

if TYPE_CHECKING:
    import numba as nb

    from fastscapelib.flow.numba.flow_kernel import (
        NumbaFlowKernel,
        NumbaFlowKernelData,
        NumbaJittedClass,
    )


class FlowKernelEroder(abc.ABC):
    """Abstract flow kernel eroder class.

    This helper class is for implementing a custom eroder based on a Numba
    :ref:`flow kernel <guide-flow-kernels>`. It has the following abstract
    methods that must be implemented in subclasses:

    - ``kernel_func`` is the kernel function
      (see also :py:func:`~fastscapelib.create_flow_kernel`)
    - ``param_spec`` and ``input_spec`` should return the names and types of the
      eroder parameters and inputs, respectively.
    - ``kernel_apply_dir`` should return the direction / order in which the kernel
      function is applied on the flow graph.

    This class provides one pre-defined "erosion" output variable for the flow
    kernel function. It also provides the same API than other eroder classes in
    Fastscapelib, i.e., an ``erode`` method that takes input variable values as
    arguments and that returns the amount erosion computed for one time step.

    Eroder parameter values should be passed as arguments of the subclass'
    ``__init__`` method, in addition to the parameters below. Those parameters
    may also be exposed as properties of the eroder subclass.

    Parameters
    ----------
    flow_graph : :py:class:`~fastscapelib.FlowGraph`
        Flow graph instance.
    **kwargs
        Keyword arguments that are passed to :py:func:`~fastscapelib.create_flow_kernel`.

    """

    _flow_graph: FlowGraph
    _kernel: NumbaFlowKernel
    _kernel_data: NumbaFlowKernelData
    _erosion: np.ndarray

    def __init__(self, flow_graph: FlowGraph, **kwargs):
        import numba as nb

        self._flow_graph = flow_graph

        kernel_sig = inspect.signature(self.kernel_func)
        if len(kernel_sig.parameters) != 1:
            raise TypeError(
                "static method 'kernel_func' must take a single argument (the kernel node data)"
            )

        spec = dict(erosion=nb.float64[::1])
        spec.update(self.input_spec())
        spec.update(self.param_spec())

        self._kernel, self._kernel_data = create_flow_kernel(
            flow_graph,
            self.kernel_func,
            spec=spec,
            outputs=["erosion"],
            apply_dir=self.kernel_apply_dir(),
            **kwargs,
        )

        self._erosion = np.zeros(flow_graph.grid_shape)
        self._kernel_data.bind(erosion=self._erosion)

    @staticmethod
    @abc.abstractmethod
    def param_spec() -> dict[str, nb.types.Type | tuple[nb.types.Type, Any]]:
        """Returns a dictionary with parameter names and their (numba) value type."""
        ...

    @staticmethod
    @abc.abstractmethod
    def input_spec() -> dict[str, nb.types.Type | tuple[nb.types.Type, Any]]:
        """Returns a dictionary with input variable names and their (numba) value type."""
        ...

    @staticmethod
    @abc.abstractmethod
    def kernel_apply_dir() -> FlowGraphTraversalDir:
        """Returns the kernel application direction and order."""
        ...

    @staticmethod
    @abc.abstractmethod
    def kernel_func(node: NumbaJittedClass) -> int | None:
        """The eroder flow kernel function."""
        ...

    def erode(self, **kwargs) -> np.ndarray:
        """Compute and returns the amount of erosion for one time step."""
        actual_inputs = set(kwargs)
        expected_inputs = set(self.input_spec())
        if missing_inputs := expected_inputs - actual_inputs:
            raise KeyError(f"inputs are missing: {missing_inputs}")
        if invalid_inputs := actual_inputs - expected_inputs:
            raise ValueError(f"invalid inputs: {invalid_inputs}")

        self._kernel_data.erosion.fill(0.0)
        self._kernel_data.bind(**kwargs)

        self._flow_graph.apply_kernel(self._kernel, self._kernel_data)

        return self._erosion.copy()

    @property
    def flow_graph(self) -> FlowGraph:
        """Returns the flow graph object used by this eroder."""
        return self._flow_graph

    @property
    def kernel(self) -> NumbaFlowKernel:
        """Returns the (Numba) flow kernel object used by this eroder."""
        return self._kernel

    @property
    def kernel_data(self) -> NumbaFlowKernelData:
        """Returns the (Numba) flow kernel data object used by this eroder."""
        return self._kernel_data
