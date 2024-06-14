import inspect

from fastscapelib.flow.numba_flow_kernel import create_flow_kernel


class NumbaEroderType(type):
    @classmethod
    def check_errs(cls, dct):
        errs = []
        missing_members = {"spec", "outputs", "apply_dir"}.difference(set(dct.keys()))
        if missing_members:
            errs.append(f"Missing mandatory class members {missing_members}")

        if "eroder_kernel" not in dct:
            errs.append(f" Missing mandatory class method 'eroder_kernel'")
        elif not isinstance(dct["eroder_kernel"], staticmethod):
            errs.append(
                "Method 'eroder_kernel' must be static (add the '@staticmethod' decorator)"
            )
        else:
            kernel_sig = inspect.signature(dct["eroder_kernel"])
            if len(kernel_sig.parameters) != 1:
                errs.append(
                    "Method 'eroder_kernel' must take a single argument (the node data)"
                )

        return errs

    def __new__(cls, name, bases, dct):
        if bases:
            errs = cls.check_errs(dct)
            if errs:
                err_msg = "\n  - ".join(errs)
                raise TypeError(f"Can't instantiate class {name}:\n  - {err_msg}")

        instance = super().__new__(cls, name, bases, dct)
        return instance


class NumbaEroderFlowKernel(metaclass=NumbaEroderType):
    def __init__(self, flow_graph, max_receivers=1, n_threads=1):
        self._flow_graph = flow_graph
        self._kernel, self._data = create_flow_kernel(
            flow_graph,
            self.eroder_kernel,
            spec=self.spec,
            outputs=self.outputs,
            max_receivers=max_receivers,
            n_threads=1,
            apply_dir=self.apply_dir,
        )

    def set_kernel_data(self, **kwargs):
        self._data.bind(**kwargs)

    @property
    def kernel_data(self):
        return self._data

    @property
    def n_threads(self):
        return self._kernel.kernel.n_threads

    @n_threads.setter
    def n_threads(self, value):
        self._kernel.kernel.n_threads = value

    @property
    def min_block_size(self):
        return self._kernel.kernel.min_block_size

    @min_block_size.setter
    def min_block_size(self, value):
        self._kernel.kernel.min_block_size = value

    @property
    def min_level_size(self):
        return self._kernel.kernel.min_level_size

    @min_level_size.setter
    def min_level_size(self, value):
        self._kernel.kernel.min_level_size = value

    def apply_kernel(self):
        self._flow_graph.apply_kernel(self._kernel, self._data)
