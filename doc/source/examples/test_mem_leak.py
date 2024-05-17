import numpy as np
import numba as nb
from testpy import Kernel, test_leak, test_no_leak


@nb.experimental.jitclass(spec=[("erosion", nb.float64[::1])])
class TestData:
    def __init__(self):
        self.erosion = np.ones(1_000_000_000)


k = Kernel()


@nb.njit
def init_node_data():
    return TestData()


jitted_fun = init_node_data.get_compile_result(
    nb.core.typing.Signature(
        TestData.class_type.instance_type,
        (),
        None,
    )
)

k.node_data_create = jitted_fun.library.get_pointer_to_function(
    jitted_fun.fndesc.llvm_cfunc_wrapper_name
)

k.node_data_free = nb.core.runtime.nrt.rtsys.library.get_pointer_to_function(
    "NRT_decref"
)

N = 4

for _ in range(N):
    test_no_leak(k)

for _ in range(N):
    test_leak(k)
