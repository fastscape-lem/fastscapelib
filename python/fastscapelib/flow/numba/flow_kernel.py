from __future__ import annotations

import ast
import inspect
import time
from collections import defaultdict
from collections.abc import Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from textwrap import dedent, indent
from types import MappingProxyType
from typing import (
    Any,
    Callable,
    Iterable,
    Iterator,
    Optional,
    ParamSpec,
    Protocol,
    Type,
    TypeVar,
    cast,
)

import numba as nb
import numpy as np
import numpy.typing as npt
from numba.experimental.jitclass import _box  # type: ignore[missing-imports]

from fastscapelib.flow import (
    FlowGraph,
    FlowGraphTraversalDir,
    _FlowKernel,
    _FlowKernelData,
)


class ConstantAssignmentVisitor(ast.NodeVisitor):
    def __init__(self, obj_name, member_name):
        self.obj_name = obj_name
        self.member_name = member_name
        self.assigned = False
        self.aliases = {obj_name}

    def visit_Assign(self, node: ast.Assign):
        """Visits a assignment nodes to search for invalid scalar variable.

        A scalar variable is passed to all node data and should not be mutated by
        a kernel.
        """
        for target in node.targets:
            if isinstance(target, ast.Attribute) and target.attr == self.member_name:
                if (
                    isinstance(target.value, ast.Name)
                    and target.value.id in self.aliases
                ):
                    self.assigned = True
            elif isinstance(target, ast.Name):
                # Track alias assignments
                if isinstance(node.value, ast.Name) and node.value.id in self.aliases:
                    self.aliases.add(target.id)
        self.generic_visit(node)

    def visit_AugAssign(self, node: ast.AugAssign):
        # Check if the augmented assignment target is an attribute of the object or its aliases
        target = node.target
        if isinstance(target, ast.Attribute) and target.attr == self.member_name:
            if isinstance(target.value, ast.Name) and target.value.id in self.aliases:
                self.assigned = True
        self.generic_visit(node)


# simplified type annotation for numba.njit decorated function
R = TypeVar("R", covariant=True)
P = ParamSpec("P")


class NumbaJittedFunc(Protocol[P, R]):
    def __call__(self, *args: P.args, **kwargs: P.kwargs) -> R: ...

    def get_compile_result(self, sig: Any) -> Any: ...


# simplified type annotation for numba.experimental.jitclass decorated class
class NumbaJittedClass(Protocol):
    class_type: Any
    _numba_type_: Any


KernelFunc = NumbaJittedFunc[[NumbaJittedClass], int]
KernelNodeDataGetter = NumbaJittedFunc[[int, NumbaJittedClass, NumbaJittedClass], int]
KernelNodeDataSetter = NumbaJittedFunc[[int, NumbaJittedClass, NumbaJittedClass], int]
KernelNodeDataCreate = NumbaJittedFunc[[], NumbaJittedClass]


@contextmanager
def timer(msg: str, do_print: bool) -> Iterator[None]:
    """Times and prints a code block."""

    start_time = time.time()
    yield
    end_time = time.time()
    elapsed_time = end_time - start_time
    if do_print:
        print(f"Time spent in {msg}: {elapsed_time:.1f} seconds")


class NumbaFlowKernelData(Mapping):
    """Proxy mapping representing the (numba) flow kernel data.

    Flow kernel data is an intermediate structure that allows accessing any
    external data (either scalar values or values defined on each node of a
    :py:class:`~fastscapelib.flow.FlowGraph`) from within the kernel function.

    It is returned by :py:func:`~fastscapelib.create_flow_kernel` and
    required by :py:meth:`~fastscapelib.FlowGraph.apply_kernel` alongside the
    :py:class:`~fastscapelib.flow.numba.flow_kernel.NumbaFlowKernel` object.

    This class implements the immutable mapping interface but still allows
    setting or updating kernel data exclusively via the ``.bind()`` method (with
    data validation).

    For convenience, this class also provides attribute-style access to kernel
    data.

    """

    _grid_shape: tuple[int, ...]
    _grid_size: np.int64
    _spec_keys: list[str]
    _grid_spec_keys: list[str]
    _bound_keys: set[str]
    _kernel_data: _FlowKernelData
    _jitclass_obj: NumbaJittedClass

    def __init__(
        self,
        grid_shape: list[int],
        spec_keys: list[str],
        grid_spec_keys: list[str],
        jitclass_obj: NumbaJittedClass,
    ):
        self._jitclass_obj = jitclass_obj
        self._grid_shape = tuple(grid_shape)
        self._grid_size = np.prod(grid_shape)
        self._spec_keys = spec_keys
        self._grid_spec_keys = grid_spec_keys
        self._bound_keys = set()

        kernel_data = _FlowKernelData()
        kernel_data.data.meminfo = _box.box_get_meminfoptr(jitclass_obj)
        kernel_data.data.data = _box.box_get_dataptr(jitclass_obj)
        self._kernel_data = kernel_data

    def __getattr__(self, name: str) -> Any:
        if name not in self._bound_keys:
            return None
        return getattr(self._jitclass_obj, name)

    def __getitem__(self, name: str) -> Any:
        if name not in self._bound_keys:
            return None
        return getattr(self._jitclass_obj, name)

    def __iter__(self):
        for key in self._spec_keys:
            yield key

    def __len__(self):
        return len(self._spec_keys)

    @property
    def jitclass_obj(self) -> NumbaJittedClass:
        """Return the numba jit-compiled class instance holding or referencing
        the kernel data.
        """
        return self._jitclass_obj

    @property
    def kernel_data(self) -> _FlowKernelData:
        """Return the object used to access kernel data from C++."""
        return self._kernel_data

    @property
    def var_names(self) -> tuple[str, ...]:
        """Return the names of all kernel input and output variables.

        Includes both scalar and array (grid) variables.

        The value of those variables may be accessed from within the kernel
        function as attributes of the node data object passed to the function as
        unique argument.

        """
        return tuple(self._spec_keys)

    @property
    def grid_var_names(self) -> tuple[str, ...]:
        """Return the names of kernel input and output array variables.

        The value of those variables may be accessed from within the kernel
        function as attributes of the node data object (passed to the function
        as unique argument) as well as attributes of the ``receivers`` and/or
        ``donors`` attributes of that node data object.

        """
        return tuple(self._grid_spec_keys)

    def bind(self, **kwargs):
        """Set or update kernel data.

        Parameters
        ----------
        **kwargs
            Argument names must be valid kernel data keys. Argument value types must correspond
            to the types defined via the kernel ``spec`` (either scalar or array).
            Kernel grid data values must be either scalar (will be expanded into an array) or
            1-d array (size matching the number of flow graph nodes) or n-d array (shape matching
            the grid shape, will be flattened).

        """
        for name, value in kwargs.items():
            if name not in self._spec_keys:
                raise KeyError(f"Unknown kernel data key {name}")

            if name in self._grid_spec_keys:
                if np.isscalar(value):
                    value = np.full(self._grid_size, value)
                if value.shape == self._grid_shape:
                    value = value.ravel()
                if value.shape != (self._grid_size,):
                    raise ValueError(
                        f"Invalid shape {value.shape} for data '{name}' (must be {(self._grid_size,)})"
                    )

            setattr(self._jitclass_obj, name, value)
            self._bound_keys.add(name)

    def check_bindings(self):
        """Check data that is bound to the kernel.

        This is called prior to executing the kernel function in order to make
        sure that all input and output data have been initialized and bound as
        kernel data.

        Raises
        ------
        ValueError
            if there's unbound data required by the kernel.

        """
        unbound_data = set(self._spec_keys) - self._bound_keys
        if unbound_data:
            raise ValueError(
                f"The following kernel data must be set prior any kernel call: {', '.join(list(unbound_data))}"
            )


@dataclass(frozen=True)
class NumbaFlowKernel:
    """Immutable Proxy object representing a numba jit-compiled flow kernel function.

    It is returned by :py:func:`~fastscapelib.create_flow_kernel` and
    required by :py:meth:`~fastscapelib.FlowGraph.apply_kernel` alongside the
    :py:class:`~fastscapelib.flow.numba.flow_kernel.NumbaFlowKernelData` mapping object.

    Attributes
    ----------
    kernel : object
        An internal object used to pass and execute the flow kernel
        either from C++ or from numba jit-compiled code.
    func : callable
        The numba-compiled kernel function, which accepts one argument
        (the kernel's node data as an instance of a numba jit-compiled class).
    apply_dir : :py:class:`~fastscapelib.FlowGraphTraversalDir`
        The direction and order in which the flow kernel will be applied along
        the graph.
    outputs : tuple
        The names of the kernel output variables.
    n_threads : int
        Number of threads to use for applying the kernel function in parallel
        along the flow graph (if equal to 1, the kernel is applied sequentially).
    generated_code : dict
        A dictionary with the source code generated for each of the
        functions below.
    generated_spec : tuple
        The numba types of the attributes of the jit-compiled classes
        used by the kernel.
    node_data_create : callable
        Internal numba compiled function used to create the kernel's
        node data object (one instance per each thread).
    node_data_init : callable or None
        Internal numba jit-compiled function used to initialize node
        data scalar values.
    node_data_getter : callable
        Internal numba jit-compiled function used to get grid data
        values at the current node, its receivers and its donors.
    node_data_setter : callable
        Internal numba jit-compiled function used to set grid data
        values from the kernel's data at the current node, its receivers
        and its donors.
    node_data_free : object
        Used internally to de-allocate memory used by the kernel's node
        data.

    """

    kernel: _FlowKernel
    func: KernelFunc
    apply_dir: FlowGraphTraversalDir
    outputs: tuple[str, ...]
    n_threads: int
    generated_code: MappingProxyType[str, str]
    generated_spec: tuple[tuple[str, Any], ...]
    node_data_create: KernelNodeDataCreate
    node_data_init: Optional[NumbaJittedFunc]
    node_data_getter: KernelNodeDataGetter
    node_data_setter: KernelNodeDataSetter
    node_data_free: Any


class NumbaFlowKernelFactory:
    """Creates a numba flow kernel.

    This factory is in charge of the heavy duty of creating a numba
    flow kernel from a kernel function and associated specs and options.

    It successively creates and compiles required jitclasses and functions
    to be used as a kernel (functions) and data.
    """

    _grid_data_ty: dict[str, nb.types.Array]
    _scalar_data_ty: dict[str, nb.types.Type]
    _generated_code: dict[str, str]
    _generated_spec: list[tuple[str, Any]]

    def __init__(
        self,
        flow_graph: FlowGraph,
        kernel_func: Callable[["NumbaJittedClass" | Any], int | None],
        spec: dict[str, nb.types.Type | tuple[nb.types.Type, Any]],
        *,
        apply_dir: FlowGraphTraversalDir = FlowGraphTraversalDir.DEPTH_UPSTREAM,
        outputs: Iterable[str] = tuple(),
        n_threads: int = 1,
        get_data_at_receivers: bool = True,
        set_data_at_receivers: bool = True,
        get_data_at_donors: bool = True,
        set_data_at_donors: bool = True,
        max_receivers: int | None = None,
        auto_resize: bool = False,
        print_stats: bool = False,
    ):
        if not outputs:
            raise ValueError("no output variable set for the flow kernel")

        if invalid_outputs := set(outputs) - set(spec):
            raise ValueError(
                f"some output variables are not defined in spec: {invalid_outputs}"
            )

        reserved_vars = {
            "receivers_count",
            "receivers_distance",
            "receivers_weight",
            "donors_count",
        }
        invalid_vars = set(spec) & reserved_vars
        if invalid_vars:
            raise ValueError(f"reserved variable names defined in spec: {invalid_vars}")

        if flow_graph.single_flow:
            max_receivers = 1
        else:
            max_receivers = flow_graph.impl().receivers.shape[1]
        max_donors = flow_graph.impl().donors.shape[1]

        with timer("flow kernel init", print_stats):
            self._flow_graph = flow_graph
            self._py_flow_kernel = kernel_func
            self._outputs = tuple(outputs)
            self._get_data_at_receivers = get_data_at_receivers
            self._set_data_at_receivers = set_data_at_receivers
            self._get_data_at_donors = get_data_at_donors
            self._set_data_at_donors = set_data_at_donors
            self._max_receivers = max_receivers
            self._max_donors = max_donors
            self._auto_resize = auto_resize
            self._spec = spec
            self._print_stats = print_stats
            self._n_threads = n_threads
            self._apply_dir = apply_dir

            self._bound_data: set[str] = set()
            self._generated_code = {}
            self._generated_spec = []

            self._set_interfaces()
            self._build_fs_kernel()

    @property
    def kernel(self):
        return NumbaFlowKernel(
            kernel=self._kernel,
            func=self.flow_kernel_func,
            apply_dir=self._apply_dir,
            outputs=self._outputs,
            n_threads=self._n_threads,
            generated_code=MappingProxyType(self._generated_code),
            generated_spec=tuple(self._generated_spec),
            node_data_create=self.node_data_create,
            node_data_init=self.node_data_init,
            node_data_getter=self.node_data_getter,
            node_data_setter=self.node_data_setter,
            node_data_free=None,
        )

    @property
    def data(self) -> NumbaFlowKernelData:
        flow_graph = self._flow_graph
        flow_graph_data = {
            "donors_idx": flow_graph.impl().donors.view(),
            "donors_count": flow_graph.impl().donors_count.view(),
            "receivers_idx": flow_graph.impl().receivers.view(),
            "receivers_count": flow_graph.impl().receivers_count.view(),
            "receivers_distance": flow_graph.impl().receivers_distance.view(),
            "receivers_weight": flow_graph.impl().receivers_weight.view(),
        }

        data = self._data_jitclass(**flow_graph_data)
        data_wrapper = NumbaFlowKernelData(
            self._flow_graph.grid_shape,
            list(self._spec.keys()),
            list(self._grid_data_ty),
            data,
        )
        data_wrapper.bind(**self._init_values)

        return data_wrapper

    def _build_fs_kernel(self):
        """Builds the fastscapelib flow kernel.

        The fastscapelib flow kernel is a collection of data and
        function pointers to be used by the `apply_kernel` runtime
        to call the flow kernel on all the grid.

        The flow kernel must be thread-safe. It can be called sequentially
        or in parallel depending on the caller implementation.
        """

        self._kernel = _FlowKernel()

        with timer("build data classes", self._print_stats):
            self._build_and_set_data_classes()
        with timer("build node data create/init/free", self._print_stats):
            self._build_and_set_node_data_create_init_free()
        with timer("build node data getter/setter", self._print_stats):
            self._build_and_set_node_data_getter()
            self._build_and_set_node_data_setter()
        with timer("build flow kernel", self._print_stats):
            self._build_and_set_flow_kernel_ptr()
        self._kernel.n_threads = self._n_threads
        self._kernel.apply_dir = self._apply_dir

    def _set_interfaces(self):
        self._grid_data_ty = {}
        self._scalar_data_ty = {}
        self._init_values = {}

        for spec_name, spec_val in self._spec.items():
            if isinstance(spec_val, tuple):
                ty, val = spec_val
            elif issubclass(spec_val.__class__, nb.types.Type):
                ty = spec_val
                val = None
            else:
                raise ValueError(
                    f"'{spec_name}' specification must be either a numba type or a (numba type, default_value) tuple"
                )

            if issubclass(ty.__class__, nb.types.Array):
                self._grid_data_ty[spec_name] = ty  # type: ignore
            else:
                self._scalar_data_ty[spec_name] = ty

            if val is not None:
                self._init_values[spec_name] = val

        invalid_outputs = [name for name in self._outputs if name not in self._spec]
        if invalid_outputs:
            raise KeyError(
                f"Output name{'s' if len(invalid_outputs)>1 else ''} {invalid_outputs} not defined in kernel specification"
            )

        invalid_outputs = [
            name for name in self._outputs if name in self._scalar_data_ty
        ]
        if invalid_outputs:
            raise TypeError(f"Output scalar data are not supported: {invalid_outputs}")

        invalid_grid_data = [
            name for name, ty in self._grid_data_ty.items() if ty.ndim != 1
        ]
        if invalid_grid_data:
            raise ValueError(
                f"Invalid shape for grid data {invalid_grid_data} (must be a flat array)"
            )

        for var in self._spec:
            if var not in self._outputs and self._is_variable_assigned_in_function(
                self._py_flow_kernel, var
            ):
                raise AttributeError(f"Invalid assignment of input variable '{var}'")

    def _build_and_set_flow_kernel_ptr(self) -> None:
        """Builds and sets the flow kernel jitted function.

        The flow kernel function is called by the computation thread with
        a node data and the integration time step.
        """
        jitted_func = cast(
            NumbaJittedFunc, nb.njit(inline="always")(self._py_flow_kernel)
        )
        self.flow_kernel_func = jitted_func

        compiled_func = jitted_func.get_compile_result(
            nb.core.typing.Signature(
                nb.none,
                (self._node_data_jitclass.class_type.instance_type,),
                None,
            )
        )
        if self._print_stats:
            print("flow kernel compilation stats", compiled_func.metadata["timers"])

        self._kernel.func = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _build_and_set_node_data_create_init_free(self):
        """Builds the node data create, init, and free functions."""

        node_data_cls = self._node_data_jitclass

        @nb.njit
        def node_data_create():
            return node_data_cls()

        self._set_node_data_create(self._node_data_jitclass, node_data_create)  # type: ignore
        self._build_and_set_node_data_init()
        self._set_node_data_free()

    def _build_and_set_node_data_init(self):
        """Sets the pointer to the node data init function.

        The node data init function is called after instantiation/creation
        to set scalar variables. Scalar variables are shared by all the nodes
        and as such are considered immutable. They can be set once before
        applying the kernel on each node.

        This function pointer may be set to a null pointer if there is no
        scalar variable to be set:
          - in case there is no scalar variable at all
          - if the scalar variable value is given in the kernel specification
            (in which case the constant scalar value is intialized inline in the
            node data constructor)
        """

        func_tmpl = dedent(
            """
        def node_data_init(node_data, data):
        {content}
        """
        )

        content = "\n".join(
            [f"node_data.{name} = data.{name}" for name in self._scalar_data_ty]
        )
        if content == "":
            self.node_data_init = None
            return

        init_source = func_tmpl.format(content=indent(content, self._indent4))

        self._generated_code["node_data_init"] = init_source

        glbls = {}
        exec(init_source, glbls)
        func = glbls["node_data_init"]

        func = cast(NumbaJittedFunc, nb.njit(inline="always", boundscheck=False)(func))
        self.node_data_init = func
        compiled_func = func.get_compile_result(
            nb.core.typing.Signature(
                nb.none,
                (
                    self._node_data_jitclass.class_type.instance_type,
                    self._data_jitclass.class_type.instance_type,
                ),
                None,
            )
        )
        if self._print_stats:
            print("Node data init compilation stats", compiled_func.metadata["timers"])

        self._kernel.node_data_init = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    @staticmethod
    def _is_variable_assigned_in_function(
        function: Callable, variable_name: str
    ) -> bool:
        source_code = inspect.getsource(function)
        tree = ast.parse(dedent(source_code))
        visitor = ConstantAssignmentVisitor(
            [*inspect.signature(function).parameters.keys()][0], variable_name
        )
        visitor.visit(tree)

        return visitor.assigned

    def _build_and_set_node_data_getter(self):
        """Builds the jitted node data getter function.

        The node data getter is called prior to the flow kernel to
        copy the global data for one specific node and its receivers
        to a node data instance.
        """

        func = self._build_node_data_getter()
        self.node_data_getter = self._build_node_data_getter()
        compiled_func = func.get_compile_result(
            nb.core.typing.Signature(
                nb.none,
                (
                    nb.uint64,
                    self._data_jitclass.class_type.instance_type,
                    self._node_data_jitclass.class_type.instance_type,
                ),
                None,
            )
        )
        if self._print_stats:
            print(
                "Node data getter compilation stats", compiled_func.metadata["timers"]
            )

        self._kernel.node_data_getter = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _build_and_set_node_data_setter(self):
        """Builds the jitted node data setter function.

        The node data setter is called after the flow kernel to
        copy back the node data in the global data instance.
        """

        self.node_data_setter = func = self._build_node_data_setter()
        compiled_func = func.get_compile_result(
            nb.core.typing.Signature(
                nb.none,
                (
                    nb.uint64,
                    self._node_data_jitclass.class_type.instance_type,
                    self._data_jitclass.class_type.instance_type,
                ),
                None,
            )
        )
        if self._print_stats:
            print(
                "Node data setter compilation stats", compiled_func.metadata["timers"]
            )

        self._kernel.node_data_setter = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _set_node_data_create(self, cls, func: NumbaJittedFunc):
        """Sets the pointer to the node data create function.

        The node data create function is called by a compute thread
        to create a new node data then used on each node of the grid.
        """

        self.node_data_create = func
        compiled_func = func.get_compile_result(
            nb.core.typing.Signature(
                cls.class_type.instance_type,
                (),
                None,
            )
        )

        self._kernel.node_data_create = compiled_func.library.get_pointer_to_function(
            compiled_func.fndesc.llvm_cfunc_wrapper_name
        )

    def _set_node_data_free(self):
        """Sets the pointer to the node data free function.

        The node data free function is called by a compute thread
        to delete a node data (it calls the jitclass destructor and
        deallocates related NRT meminfo and data pointers).
        """

        self._kernel.node_data_free = (
            nb.core.runtime.nrt.rtsys.library.get_pointer_to_function("NRT_decref")
        )

    def _build_and_set_data_classes(self):
        """Builds data and node data jitclasses."""

        self._build_node_data_jitclass()
        self._build_data_jitclass()

    @staticmethod
    def _generate_jitclass(
        name,
        spec: Iterable[tuple[str, nb.types.Type]],
        init_source: str,
        glbls: dict = {},
    ) -> Type[NumbaJittedClass]:
        exec(init_source, glbls)
        ctor = glbls["generated_init"]
        return cast(
            Type[NumbaJittedClass],
            nb.experimental.jitclass(spec)(type(name, (), {"__init__": ctor})),
        )

    def _build_node_data_jitclass(self):
        """Builds a node data jitclass.

        A node data instance contains all the data related to a node and its receivers.
        It is used to call the flow kernel function with this specific node data.

        The grid data specified by the user are aggregated with the constants and
        also with the same data at the receivers. The count, distance and weight to
        the current node's receivers are also provided.

        We need to workaround some `numba` limitations:
        - `setattr` is not implemented -> use a template source code to be executed
        """

        grid_data = self._grid_data_ty
        scalar_data = self._scalar_data_ty
        grid_data_ty_spec = [(name, ty) for name, ty in grid_data.items()]

        receivers_flow_graph_spec = [
            ("distance", nb.float64[::1]),
            ("weight", nb.float64[::1]),
        ]
        receivers_internal_spec = [
            ("count", nb.uint64),
        ]

        receivers_spec = (
            receivers_flow_graph_spec
            + grid_data_ty_spec
            + receivers_internal_spec
            + [
                ("_" + name, ty)
                for name, ty in receivers_flow_graph_spec + grid_data_ty_spec
            ]
        )

        @nb.experimental.jitclass(receivers_spec)
        class ReceiversData(object):
            def __init__(self):
                pass

        donors_internal_spec = [
            ("count", nb.uint64),
        ]
        donors_spec = (
            grid_data_ty_spec
            + donors_internal_spec
            + [("_" + name, ty) for name, ty in grid_data_ty_spec]
        )

        @nb.experimental.jitclass(donors_spec)
        class DonorsData(object):
            def __init__(self):
                pass

        base_spec = [
            ("receivers", ReceiversData.class_type.instance_type),
            ("donors", DonorsData.class_type.instance_type),
        ]
        grid_data_spec = [(name, ty.dtype) for name, ty in grid_data.items()]
        scalar_data_spec = [(name, ty) for name, ty in scalar_data.items()]

        __init___template = dedent(
            """
        def generated_init(self):
            self.receivers = ReceiversData()
            self.receivers._distance = np.ones({default_rec_size})
            self.receivers._weight = np.ones({default_rec_size})
            {receivers_content_init}

            self.receivers.distance = self.receivers._distance[:]
            self.receivers.weight = self.receivers._weight[:]
            {receivers_content_view}

            self.receivers.count = 0

            self.donors = DonorsData()
            {donors_content_init}

            {donors_content_view}

            self.donors.count = 0
        """
        )

        default_rec_size = 0 if self._auto_resize else self._max_receivers
        receivers_content_init = "\n    ".join(
            [
                f"self.receivers._{name} = np.ones({default_rec_size}, dtype=np.{value.dtype})"
                for name, value in grid_data.items()
            ]
        )
        receivers_content_view = "\n    ".join(
            [f"self.receivers.{name} = self.receivers._{name}[:]" for name in grid_data]
        )

        default_don_size = 0 if self._auto_resize else self._max_donors
        donors_content_init = "\n    ".join(
            [
                f"self.donors._{name} = np.ones({default_don_size}, dtype=np.{value.dtype})"
                for name, value in grid_data.items()
            ]
        )
        donors_content_view = "\n    ".join(
            [f"self.donors.{name} = self.donors._{name}[:]" for name in grid_data]
        )

        init_source = __init___template.format(
            default_rec_size=default_rec_size,
            receivers_content_init=receivers_content_init,
            receivers_content_view=receivers_content_view,
            default_don_size=default_don_size,
            donors_content_init=donors_content_init,
            donors_content_view=donors_content_view,
        )

        self._generated_spec = base_spec + grid_data_spec + scalar_data_spec
        self._generated_code["node_data_jitclass_init"] = init_source

        self._node_data_jitclass = NumbaFlowKernelFactory._generate_jitclass(
            "FlowKernelNodeData",
            self._generated_spec,
            init_source,
            {"ReceiversData": ReceiversData, "DonorsData": DonorsData, "np": np},
        )

    def _build_data_jitclass(self):
        """Builds a data jitclass.

        A data instance contains all the data maped on the grid. It contains the
        initial and final data after applying the kernel.

        The grid data specified by the user are aggregated with the flow graph ones
        to give access to the receivers (indices, count, distance, weight).
        The flow graph donors are also exposed but currently not used.
        """

        grid_data_ty = self._grid_data_ty
        scalar_data_ty = self._scalar_data_ty

        base_spec = dict(
            donors_idx=nb.uint64[:, ::1],
            donors_count=nb.uint64[::1],
            receivers_idx=nb.uint64[:, ::1],
            receivers_count=nb.uint64[::1],
            receivers_distance=nb.float64[:, ::1],
            receivers_weight=nb.float64[:, ::1],
        )
        grid_data_spec = [(name, ty) for name, ty in grid_data_ty.items()]
        scalar_data_spec = [(name, ty) for name, ty in scalar_data_ty.items()]
        spec = (
            [(name, ty) for name, ty in base_spec.items()]
            + grid_data_spec
            + scalar_data_spec
        )

        __init___template = dedent(
            """
        def generated_init(self, {args}):
            {content}
        """
        )

        content = "\n    ".join([f"self.{name} = {name}" for name in base_spec])
        args = ", ".join([name for name in base_spec])
        init_source = __init___template.format(content=content, args=args)

        self._generated_code["data_jitclass_init"] = init_source

        self._data_jitclass = NumbaFlowKernelFactory._generate_jitclass(
            "FlowKernelData",
            spec,
            init_source,
        )

    _node_data_getter_tmpl = dedent(
        """
        def node_data_getter(index, data, node_data):
            # --- data at node ---
        {node_content}

            # --- data at flow receivers ---
            receivers_count = data.receivers_count[index]
            receivers = node_data.receivers

        {receivers_resize_content}

            receivers.count = receivers_count

            for i in range(receivers_count):
                receiver_idx = data.receivers_idx[index, i]
                receivers._distance[i] = data.receivers_distance[index, i]
                receivers._weight[i] = data.receivers_weight[index, i]
        {receivers_set_content}

            # --- data at flow donors ---
            donors_count = data.donors_count[index]
            donors = node_data.donors

        {donors_resize_content}

            donors.count = donors_count

            for i in range(donors_count):
                donor_idx = data.donors_idx[index, i]
        {donors_set_content}

            return 0
        """
    )

    _node_data_getter_fixed_resize_tmpl = dedent(
        """\
        if {receivers_or_donors}_count > {max_count}:
            return 1

        if {receivers_or_donors}_count != {receivers_or_donors}.count:
        {set_views}
        """
    ).rstrip("\n")

    _node_data_getter_dynamic_resize_tmpl = dedent(
        """\
        if {receivers_or_donors}_count != {receivers_or_donors}.count:
            if {receivers_or_donors}_count > {receivers_or_donors}.count:
        {resize_source}
        {set_views}
        """
    ).rstrip("\n")

    _set_view_tmpl = dedent(
        """
        set_view(
            (
        {view_data},
            ),
            {receivers_or_donors}_count
        )"""
    ).lstrip("\n")

    _indent4 = " " * 4
    _indent8 = _indent4 * 2
    _indent12 = _indent4 * 3

    def _build_node_data_getter(self) -> NumbaJittedFunc:
        """Builds a node data getter from the global data

        The node data getter is called prior to the flow kernel to
        copy the global data for one specific node and its receivers
        to a node data instance.

        We need to workaround some `numba` limitations:
        - `setattr` is not implemented -> use a template source code to be executed
        - passing a tuple of various array types doesn't work -> use `set_view` with
          pack of consistent data

        Note: instead of using `set_view` multiple times, an inlining of the views has
        been tested but generates very poor performances (not understood in details)
        """

        if self._auto_resize:
            resize_tmpl = self._node_data_getter_dynamic_resize_tmpl
        else:
            resize_tmpl = self._node_data_getter_fixed_resize_tmpl

        # --- node data ---
        node_content = "\n".join(
            [f"node_data.{name} = data.{name}[index]" for name in self._grid_data_ty]
        )

        # --- flow receivers data ---
        receivers_grid_data_ty = {
            "distance": nb.float64[::1],
            "weight": nb.float64[::1],
        }
        if self._get_data_at_receivers:
            receivers_grid_data_ty.update(self._grid_data_ty.copy())

        receivers_data_dtypes: dict[nb.types.Type, list[str]] = defaultdict(list)
        for name, value in receivers_grid_data_ty.items():
            receivers_data_dtypes[value.dtype].append(name)

        receivers_view_data = [
            f",\n".join([f"(receivers.{name}, receivers._{name})" for name in names])
            for names in receivers_data_dtypes.values()
        ]

        receivers_resize_source = "\n".join(
            [
                f"receivers._{name} = np.empty(receivers_count, dtype=np.{value.dtype})"
                for name, value in receivers_grid_data_ty.items()
            ]
        )

        if self._get_data_at_receivers:
            receivers_set_content = "\n".join(
                [
                    f"receivers._{name}[i] = data.{name}[receiver_idx]"
                    for name in self._grid_data_ty
                ]
            )
        else:
            receivers_set_content = ""

        receivers_resize_content = resize_tmpl.format(
            set_views=indent(
                "\n".join(
                    self._set_view_tmpl.format(
                        view_data=indent(content, self._indent8),
                        receivers_or_donors="receivers",
                    )
                    for content in receivers_view_data
                ),
                self._indent4,
            ),
            resize_source=indent(receivers_resize_source, self._indent8),
            receivers_or_donors="receivers",
            max_count=self._max_receivers,
        )

        # --- flow donors data ---
        if not self._get_data_at_donors:
            donors_resize_content = ""
            donors_set_content = ""
        else:
            donors_grid_data_ty = self._grid_data_ty.copy()

            donors_data_dtypes: dict[nb.types.Type, list[str]] = defaultdict(list)
            for name, value in donors_grid_data_ty.items():
                donors_data_dtypes[value.dtype].append(name)

            donors_view_data = [
                f",\n".join([f"(donors.{name}, donors._{name})" for name in names])
                for names in donors_data_dtypes.values()
            ]

            donors_resize_source = "\n".join(
                [
                    f"donors._{name} = np.empty(donors_count, dtype=np.{value.dtype})"
                    for name, value in donors_grid_data_ty.items()
                ]
            )

            donors_set_content = "\n".join(
                [
                    f"donors._{name}[i] = data.{name}[donor_idx]"
                    for name in self._grid_data_ty
                ]
            )

            donors_resize_content = resize_tmpl.format(
                set_views=indent(
                    "\n".join(
                        self._set_view_tmpl.format(
                            view_data=indent(content, self._indent8),
                            receivers_or_donors="donors",
                        )
                        for content in donors_view_data
                    ),
                    self._indent4,
                ),
                resize_source=indent(donors_resize_source, self._indent8),
                receivers_or_donors="donors",
                max_count=self._max_donors,
            )

        getter_source = self._node_data_getter_tmpl.format(
            node_content=indent(node_content, self._indent4),
            receivers_resize_content=indent(receivers_resize_content, self._indent4),
            receivers_set_content=indent(receivers_set_content, self._indent8),
            donors_resize_content=indent(donors_resize_content, self._indent4),
            donors_set_content=indent(donors_set_content, self._indent8),
        )
        self._generated_code["node_data_getter"] = getter_source

        @nb.njit(inline="always")
        def set_view(data, size):
            for view, source in data:
                view = source[:size]

        glbls = {"np": np, "set_view": set_view}
        exec(getter_source, glbls)
        getter_fun = glbls["node_data_getter"]

        return nb.njit(inline="always", boundscheck=False)(getter_fun)

    node_data_setter_tmpl = dedent(
        """
        def node_data_setter(index, node_data, data):
            # --- data at flow receivers ---
            receivers_count = data.receivers_count[index]
            receivers = node_data.receivers

            for i in range(receivers_count):
                receiver_idx = data.receivers_idx[index, i]
        {receivers_content}

            # --- data at flow donors ---
            donors_count = data.donors_count[index]
            donors = node_data.donors

            for i in range(donors_count):
                donor_idx = data.donors_idx[index, i]
        {donors_content}

            # --- data at node ---
            # Note: at a base level node (receiver index == node index):
            # if a new value was set above via the receiver, it will
            # be overwritten here.
        {node_content}

            return 0
        """
    )

    def _build_node_data_setter(self) -> NumbaJittedFunc:
        """Builds a node data setter from the global data

        The node data setter is called after the flow kernel to
        copy back the node data in the global data instance.

        We need to workaround some `numba` limitations:
        - `setattr` is not implemented -> use a template source code to be executed
        """

        node_content = "\n".join(
            [
                f"data.{name}[index] = node_data.{name}"
                for name in self._grid_data_ty
                if name in self._outputs
            ]
        )

        if self._set_data_at_receivers:
            receivers_content = "\n".join(
                [
                    f"data.{name}[receiver_idx] = receivers.{name}[i]"
                    for name in self._grid_data_ty
                    if name in self._outputs
                ]
            )
        else:
            receivers_content = ""

        if self._set_data_at_donors:
            donors_content = "\n".join(
                [
                    f"data.{name}[donor_idx] = donors.{name}[i]"
                    for name in self._grid_data_ty
                    if name in self._outputs
                ]
            )
        else:
            donors_content = ""

        setter_source = self.node_data_setter_tmpl.format(
            node_content=indent(node_content, self._indent4),
            receivers_content=indent(receivers_content, self._indent8),
            donors_content=indent(donors_content, self._indent8),
        )

        self._generated_code["node_data_setter"] = setter_source

        glbls: dict[str, Any] = dict()
        exec(setter_source, glbls)
        setter_fun = glbls["node_data_setter"]
        return nb.njit(inline="always", boundscheck=False)(setter_fun)


@nb.njit
def _apply_flow_kernel(
    indices: np.ndarray,
    func: KernelFunc,
    data: NumbaJittedClass,
    node_data: NumbaJittedClass,
    node_data_getter: KernelNodeDataGetter,
    node_data_setter: KernelNodeDataSetter,
):
    """Applies a kernel on a grid.

    This Python implementation allows numba to inline
    the calls to the node data getter, the flow kernel,
    and node data setter.

    In case of sequential execution, it may give better
    performances than using the C++ fastscapelib `apply_kernel`
    implementation.
    """

    for i in indices:
        if node_data_getter(i, data, node_data):
            raise ValueError(
                f"Invalid index {i} encountered in node_data getter function. "
                "Please check if you are using dynamic receivers count "
                "('max_receivers=None') or adjust this setting in "
                "`create_flow_kernel()`."
            )
        func(node_data)
        node_data_setter(i, node_data, data)

    return 0


def apply_flow_kernel(
    flow_graph: FlowGraph, kernel: NumbaFlowKernel, data: NumbaFlowKernelData
):
    """Applies a kernel on a grid.

    This wrapper function calls the Python jitted implementation
    of a sequential call of the flow kernel on the grid nodes.
    """

    wrapped_kernel = kernel.kernel
    node_data = kernel.node_data_create()
    if kernel.node_data_init:
        kernel.node_data_init(node_data, data.jitclass_obj)

    indices: npt.NDArray[np.uint64]

    if wrapped_kernel.apply_dir == FlowGraphTraversalDir.ANY:
        indices = np.arange(0, flow_graph.size, 1, dtype=np.uint64)
    elif wrapped_kernel.apply_dir == FlowGraphTraversalDir.BREADTH_UPSTREAM:
        indices = flow_graph.impl().bfs_indices
    elif wrapped_kernel.apply_dir == FlowGraphTraversalDir.DEPTH_UPSTREAM:
        indices = flow_graph.impl().dfs_indices
    elif wrapped_kernel.apply_dir == FlowGraphTraversalDir.BREADTH_DOWNSTREAM:
        indices = flow_graph.impl().bfs_indices[::-1]
    elif wrapped_kernel.apply_dir == FlowGraphTraversalDir.DEPTH_DOWNSTREAM:
        indices = flow_graph.impl().dfs_indices[::-1]
    else:
        raise ValueError(
            f"Unknown kernel application direction: {wrapped_kernel.apply_dir!r}"
        )

    return _apply_flow_kernel(
        indices,
        kernel.func,
        data.jitclass_obj,
        node_data,
        kernel.node_data_getter,
        kernel.node_data_setter,
    )
