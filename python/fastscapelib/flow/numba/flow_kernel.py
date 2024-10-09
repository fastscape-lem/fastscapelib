from __future__ import annotations

import ast
import inspect
import time
from collections import defaultdict
from collections.abc import Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from textwrap import dedent, indent
from typing import (
    TYPE_CHECKING,
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
    def __call__(self, *args: P.args, **kwargs: P.kwargs) -> R:
        ...

    def get_compile_result(self, sig: Any) -> Any:
        ...


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
    """Proxy mapping for access to numba flow kernel data.

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
        unbound_data = set(self._spec_keys) - self._bound_keys
        if unbound_data:
            raise ValueError(
                f"The following kernel data must be set prior any kernel call: {', '.join(list(unbound_data))}"
            )


@dataclass
class NumbaFlowKernel:
    """Stores a numba flow kernel."""

    kernel: _FlowKernel
    node_data_create: KernelNodeDataCreate
    node_data_init: Optional[NumbaJittedFunc]
    node_data_getter: KernelNodeDataGetter
    node_data_setter: KernelNodeDataSetter
    node_data_free: Any
    func: KernelFunc


class NumbaFlowKernelFactory:
    """Creates a numba flow kernel.

    This factory is in charge of the heavy duty of creating a numba
    flow kernel from a kernel function and associated specs and options.

    It successively creates and compiles required jitclasses and functions
    to be used as a kernel (functions) and data.
    """

    _grid_data_ty: dict[str, nb.types.Array]
    _scalar_data_ty: dict[str, nb.types.Type]

    def __init__(
        self,
        flow_graph: FlowGraph,
        kernel_func: Callable[["NumbaJittedClass"], int],
        spec: dict[str, nb.types.Type | tuple[nb.types.Type, Any]],
        apply_dir: FlowGraphTraversalDir,
        outputs: Iterable[str] = tuple(),
        max_receivers: int = -1,
        n_threads: int = 1,
        print_generated_code: bool = False,
        print_stats: bool = False,
    ):
        with timer("flow kernel init", print_stats):
            self._flow_graph = flow_graph
            self._py_flow_kernel = kernel_func
            self._outputs = outputs
            self._max_receivers = max_receivers
            self._spec = spec
            self._print_generated_code = print_generated_code
            self._print_stats = print_stats
            self._n_threads = n_threads
            self._apply_dir = apply_dir

            self._bound_data: set[str] = set()

            self._set_interfaces()
            self._build_fs_kernel()

    @property
    def kernel(self):
        return NumbaFlowKernel(
            kernel=self._kernel,
            node_data_create=self.node_data_create,
            node_data_init=self.node_data_init,
            node_data_getter=self.node_data_getter,
            node_data_setter=self.node_data_setter,
            node_data_free=None,
            func=self.flow_kernel_func,
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

        if self._print_generated_code:
            print(f"Node data init source code:\n{indent(init_source, self._indent4)}")

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

        func = self._build_node_data_getter(self._max_receivers)
        self.node_data_getter = self._build_node_data_getter(self._max_receivers)
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
        max_receivers = self._max_receivers

        receivers_flow_graph_spec = [
            ("distance", nb.float64[::1]),
            ("weight", nb.float64[::1]),
        ]
        receivers_grid_data_ty_spec = [(name, ty) for name, ty in grid_data.items()]
        receivers_internal_spec = [
            ("count", nb.uint64),
        ]

        receivers_spec = (
            receivers_flow_graph_spec
            + receivers_grid_data_ty_spec
            + receivers_internal_spec
            + [
                ("_" + name, ty)
                for name, ty in receivers_flow_graph_spec + receivers_grid_data_ty_spec
            ]
        )

        @nb.experimental.jitclass(receivers_spec)
        class ReceiversData(object):
            def __init__(self):
                pass

        base_spec = [("receivers", ReceiversData.class_type.instance_type)]
        grid_data_spec = [(name, ty.dtype) for name, ty in grid_data.items()]
        scalar_data_spec = [(name, ty) for name, ty in scalar_data.items()]

        __init___template = dedent(
            """
        def generated_init(self):
            self.receivers = ReceiversData()
            self.receivers._distance = np.ones({default_size})
            self.receivers._weight = np.ones({default_size})
            {receivers_content_init}

            self.receivers.distance = self.receivers._distance[:]
            self.receivers.weight = self.receivers._weight[:]
            {receivers_content_view}

            self.receivers.count = 0
        """
        )

        default_size = max_receivers if max_receivers > 0 else 0
        receivers_content_init = "\n    ".join(
            [
                f"self.receivers._{name} = np.ones({default_size}, dtype=np.{value.dtype})"
                for name, value in grid_data.items()
            ]
        )
        receivers_content_view = "\n    ".join(
            [f"self.receivers.{name} = self.receivers._{name}[:]" for name in grid_data]
        )
        init_source = __init___template.format(
            receivers_content_init=receivers_content_init,
            receivers_content_view=receivers_content_view,
            default_size=default_size,
        )

        spec = base_spec + grid_data_spec + scalar_data_spec

        if self._print_generated_code:
            print(f"Node data jitclass spec:\n{spec}")
            print(
                f"Node data jitclass constructor source code:\n{indent(init_source, self._indent4)}"
            )

        self._node_data_jitclass = NumbaFlowKernelFactory._generate_jitclass(
            "FlowKernelNodeData",
            base_spec + grid_data_spec + scalar_data_spec,
            init_source,
            {"ReceiversData": ReceiversData, "np": np},
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

        if self._print_generated_code:
            print(
                f"Data jitclass constructor source code:\n{indent(init_source, self._indent4)}"
            )

        self._data_jitclass = NumbaFlowKernelFactory._generate_jitclass(
            "FlowKernelData",
            spec,
            init_source,
        )

    _node_data_getter_tmpl = dedent(
        """
        def node_data_getter(index, data, node_data):
            receivers_count = data.receivers_count[index]
            receivers = node_data.receivers

        {node_content}
            {resize_content}

            receivers.count = receivers_count

            for i in range(receivers_count):
                receiver_idx = data.receivers_idx[index, i]
                receivers._distance[i] = data.receivers_distance[index, i]
                receivers._weight[i] = data.receivers_weight[index, i]
        {receivers_set_content}

            return 0
        """
    )

    _node_data_getter_fixed_resize_tmpl = dedent(
        """
        if {max_receivers} < receivers_count:
            return 1

        if receivers_count != receivers.count:
        {set_views}
        """
    ).rstrip("\n")

    _node_data_getter_dynamic_resize_tmpl = dedent(
        """
        if receivers_count != receivers.count:
            if receivers_count > receivers.count:
        {receivers_resize_source}
        {set_views}
        """
    ).rstrip("\n")

    _set_view_tmpl = dedent(
        """
        set_view(
            (
        {view_data},
            ),
            receivers_count
        )"""
    ).lstrip("\n")

    _indent4 = " " * 4
    _indent8 = _indent4 * 2
    _indent12 = _indent4 * 3

    def _build_node_data_getter(self, max_receivers: int) -> NumbaJittedFunc:
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

        receivers_grid_data_ty = self._grid_data_ty.copy()
        receivers_grid_data_ty.update(
            {"distance": nb.float64[::1], "weight": nb.float64[::1]}
        )

        data_dtypes: dict[nb.types.Type, list[str]] = defaultdict(list)
        for name, value in receivers_grid_data_ty.items():
            data_dtypes[value.dtype].append(name)

        node_content = "\n".join(
            [f"node_data.{name} = data.{name}[index]" for name in self._grid_data_ty]
        )

        receivers_view_data = [
            f",\n".join([f"(receivers.{name}, receivers._{name})" for name in names])
            for names in data_dtypes.values()
        ]

        receivers_resize_source = "\n".join(
            [
                f"receivers._{name} = np.empty(receivers_count, dtype=np.{value.dtype})"
                for name, value in receivers_grid_data_ty.items()
            ]
        )
        receivers_set_content = "\n".join(
            [
                f"receivers._{name}[i] = data.{name}[receiver_idx]"
                for name in self._grid_data_ty
            ]
        )

        if max_receivers > 0:
            resize_tmpl = self._node_data_getter_fixed_resize_tmpl
        else:
            resize_tmpl = self._node_data_getter_dynamic_resize_tmpl

        resize_source = indent(
            resize_tmpl.format(
                set_views=indent(
                    "\n".join(
                        self._set_view_tmpl.format(
                            view_data=indent(content, self._indent8)
                        )
                        for content in receivers_view_data
                    ),
                    self._indent4,
                ),
                receivers_resize_source=indent(receivers_resize_source, self._indent8),
                max_receivers=max_receivers,
            ),
            self._indent4,
        )

        getter_source = self._node_data_getter_tmpl.format(
            node_content=indent(node_content, self._indent4),
            resize_content=resize_source,
            receivers_set_content=indent(receivers_set_content, self._indent8),
        )
        if self._print_generated_code:
            print(
                f"Node data getter source code:\n{indent(getter_source, self._indent4)}"
            )

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
        {content}

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

        content = "\n".join(
            [
                f"data.{name}[index] = node_data.{name}"
                for name in self._grid_data_ty
                if name in self._outputs
            ]
        )

        setter_source = self.node_data_setter_tmpl.format(
            content=indent(content, self._indent4)
        )
        if self._print_generated_code:
            print(
                f"Node data setter source code:\n{indent(setter_source, self._indent4)}"
            )

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
                f"Invalid index {i} encountered in node_data getter function\n"
                "Please check if you are using dynamic receivers count "
                "('max_receivers=-1') or adjust this setting in the "
                "'Kernel' specification"
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

    if wrapped_kernel.apply_dir == FlowGraphTraversalDir.ANY:
        indices = np.arange(0, flow_graph.size, 1, dtype=np.uint64)
    elif wrapped_kernel.apply_dir == FlowGraphTraversalDir.BREADTH_UPSTREAM:
        indices = flow_graph.impl().bfs_indices
    elif wrapped_kernel.apply_dir == FlowGraphTraversalDir.DEPTH_UPSTREAM:
        indices = flow_graph.impl().dfs_indices
    else:
        raise ValueError("Unsupported kernel application direction")

    return _apply_flow_kernel(
        indices,
        kernel.func,
        data.jitclass_obj,
        node_data,
        kernel.node_data_getter,
        kernel.node_data_setter,
    )
