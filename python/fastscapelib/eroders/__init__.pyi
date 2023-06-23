from typing import overload

import numpy as np
import numpy.typing as npt

from fastscapelib.flow import FlowGraph
from fastscapelib.grid import RasterGrid

class SPLEroder:
    @overload
    def __init__(
        self,
        flow_graph: FlowGraph,
        k_coef: float,
        area_exp: float,
        slope_exp: float,
        tolerance: float,
    ) -> None: ...
    @overload
    def __init__(
        self,
        flow_graph: FlowGraph,
        k_coef: npt.NDArray[np.float64],
        area_exp: float,
        slope_exp: float,
        tolerance: float,
    ) -> None: ...
    @property
    def k_coef(self) -> npt.ArrayLike: ...
    @k_coef.setter
    def k_coef(self, value: npt.ArrayLike) -> None: ...
    @property
    def area_exp(self) -> float: ...
    @area_exp.setter
    def area_exp(self, value: float) -> None: ...
    @property
    def slope_exp(self) -> float: ...
    @slope_exp.setter
    def slope_exp(self, value: float) -> None: ...
    @property
    def tolerance(self) -> float: ...
    @property
    def n_corr(self) -> int: ...
    def erode(
        self,
        elevation: npt.NDArray[np.float64],
        drainage_area: npt.NDArray[np.float64],
        dt: float,
    ) -> npt.NDArray[np.float64]: ...

class DiffusionADIEroder:
    @overload
    def __init__(self, grid: RasterGrid, k_coef: float) -> None: ...
    @overload
    def __init__(self, grid: RasterGrid, k_coef: npt.NDArray[np.float64]) -> None: ...
    @property
    def k_coef(self) -> npt.ArrayLike: ...
    @k_coef.setter
    def k_coef(self, value: npt.ArrayLike) -> None: ...
    def erode(
        self, elevation: npt.NDArray[np.float64], dt: float
    ) -> npt.NDArray[np.float64]: ...
