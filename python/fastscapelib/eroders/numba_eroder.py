import inspect
from fastscapelib.flow import create_flow_kernel


class NumbaEroderType(type):

    @classmethod
    def check_errs(cls, dct):
        errs = []
        missing_members = {"spec", "outputs", "application_order"}.difference(set(dct.keys()))
        if missing_members:
            errs.append(f"Missing mandatory class members {missing_members}")
        
        if "eroder_kernel" not in dct:
            errs.append(f" Missing mandatory class method 'eroder_kernel'")
        elif not isinstance(dct["eroder_kernel"], staticmethod):
            errs.append("Method 'eroder_kernel' must be static (add the '@staticmethod' decorator)")
        else:
            kernel_sig = inspect.signature(dct["eroder_kernel"])
            if len(kernel_sig.parameters) != 1:
                errs.append("Method 'eroder_kernel' must take a single argument (the node data)")

        return errs
    
    def __new__(cls, name, bases, dct):
        if bases:
            errs = cls.check_errs(dct)
            if errs:
                raise TypeError(f"Can't instantiate class {name}:\n  - {'\n  - '.join(errs)}")  
    
        instance = super().__new__(cls, name, bases, dct)
        return instance
        

class NumbaEroderBase(metaclass=NumbaEroderType):

    def __init__(self, flow_graph, max_receivers=1, n_threads=1):
        self._flow_graph = flow_graph
        self._kernel, self._data = create_flow_kernel(
            flow_graph,
            self.eroder_kernel,
            spec=self.spec,
            outputs=self.outputs,
            max_receivers=max_receivers,
            n_threads=1,
            application_order=self.application_order
        )

    def set_data(self, **kwargs):
        self._data.bind(**kwargs)

    @property
    def data(self):
        return self._data

    @property
    def n_threads(self):
        return self._kernel.kernel.n_threads

    @n_threads.setter
    def n_threads(self, value):
        self._kernel.kernel.n_threads = value

    def apply_kernel(self):
        self._flow_graph.apply_kernel(self._kernel, self._data)
