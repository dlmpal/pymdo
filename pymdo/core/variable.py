from os.path import join
from typing import Dict
from typing import List
from typing import Callable
from dataclasses import dataclass
from dataclasses import field

from numpy import float64
from numpy import complex128
from numpy import inf
from numpy import zeros
from numpy import ndarray
from numpy import mean
from numpy import isinf
from numpy import isneginf
from numpy.linalg import norm
import matplotlib.pyplot as plt

FLOAT_DATA_TYPE = float64
""" Datatype for floating point arithmetic """

COMPLEX_FLOAT_DATA_TYPE = complex128
""" Datatype for complex floating point arithmetic. 
    Needs to be twice the size of FLOAT_DATA_TYPE 
"""


@dataclass(frozen=True)
class Variable:

    name: str = field(default=None, hash=True)

    size: int = field(default=1, hash=False)

    lb: float = field(default=-inf, hash=False)

    ub: float = field(default=inf, hash=False)

    keepFeasible: bool = field(default=True, hash=False)

    def norm_values(self, _val: ndarray) -> ndarray:
        return (_val - self.lb) / (self.ub - self.lb)

    def norm_grad(self, _grad: ndarray) -> ndarray:
        return _grad * (self.ub - self.lb)

    def denorm_values(self, _val: ndarray) -> ndarray:
        return _val * (self.ub - self.lb) + self.lb

    def denorm_grad(self, _grad: ndarray) -> ndarray:
        return _grad / (self.ub - self.lb)


def dict_to_array(vars: List[Variable], valueDict: Dict[str, ndarray]) -> ndarray:

    nVars = sum([var.size for var in vars])

    valueArray = zeros(nVars, FLOAT_DATA_TYPE)

    r = 0

    for var in vars:

        valueArray[r: r + var.size] = valueDict[var.name]

        r += var.size

    return valueArray


def array_to_dict(vars: List[Variable], valueArray: ndarray) -> Dict[str, ndarray]:

    valueDict = {}

    r = 0

    for var in vars:

        valueDict[var.name] = valueArray[r: r + var.size]

        r += var.size

    return valueDict


def dict_to_array2d(inputVars: List[Variable], outputVars: List[Variable], valueDict: Dict[str, Dict[str, ndarray]],
                    flatten: bool = False) -> ndarray:

    nInVars = sum([var.size for var in inputVars])

    nOutVars = sum([var.size for var in outputVars])

    valueArray = zeros((nOutVars, nInVars),
                       FLOAT_DATA_TYPE)

    r = 0

    for Fi in outputVars:

        c = 0

        for xj in inputVars:

            valueArray[r: r + Fi.size,
                       c: c + xj.size] = valueDict[Fi.name][xj.name]

            c += xj.size

        r += Fi.size

    if flatten:
        return valueArray.reshape(-1)
    return valueArray


def array2d_to_dict(inputVars: List[Variable], outputVars: List[Variable], valueArray: ndarray) -> Dict[str, Dict[str, ndarray]]:

    varDict: Dict[str, Dict[str, ndarray]] = {}

    r = 0

    for Fi in outputVars:

        c = 0

        varDict[Fi.name]: Dict[str, ndarray] = {}

        for xj in inputVars:

            varDict[Fi.name][xj.name] = valueArray[r: r +
                                                   Fi.size, c: c + xj.size]

            c += xj.size

        r += Fi.size

    return varDict


def normalize_design_vector(vars: List[Variable], deisgnVector: Dict[str, ndarray]) -> Dict[str, ndarray]:

    for var in vars:

        deisgnVector[var.name] = var.norm_values(deisgnVector[var.name])

    return deisgnVector


def denormalize_design_vector(vars: List[Variable], designVector: Dict[str, ndarray]) -> Dict[str, ndarray]:

    for var in vars:

        designVector[var.name] = var.denorm_values(designVector[var.name])

    return designVector


def normalize_gradient(inputVars: List[Variable], outputVars: List[Variable], grad: Dict[str, Dict[str, ndarray]]):

    for outVar in outputVars:

        for inVar in inputVars:

            grad[outVar.name][inVar.name] = inVar.norm_grad(
                grad[outVar.name][inVar.name])

    return grad


def denormalize_gradient(inputVars: List[Variable], outputVars: List[Variable], grad: Dict[str, Dict[str, ndarray]]):

    for outVar in outputVars:

        for inVar in inputVars:

            grad[outVar.name][inVar.name] = inVar.denorm_grad(
                grad[outVar.name][inVar.name])

    return grad


def check_values_match(vars: List[Variable],
                       original: Dict[str, ndarray] = {},
                       test: Dict[str, ndarray] = {},
                       tol: float = 1e-9) -> bool:

    if not original or not test:

        return False

    err = 0.0

    for var in vars:

        err += (norm(original[var.name] - test[var.name]
                     ) / (1.0 + norm(test[var.name])))

    if err >= tol:

        return False

    return True


def single_variable_plot(
    var: Variable,
    values: List[ndarray],
    with_bounds: bool = True,
    with_legend: bool = True,
    with_grid: bool = True,
    with_labels: bool = True,
    xlabel: str = "iterations",
    ylabel: str = None,
    show: bool = True,
    save: bool = False,
    path: str = ""
) -> None:

    fig, ax = plt.subplots(nrows=1,
                           ncols=1)

    xValues = [i for i in range(len(values))]

    yValues = [[values[j][i]
                for j in range(len(values))] for i in range(var.size)]

    for i in range(var.size):

        plt.plot(xValues,
                 yValues[i],
                 label=f"{var.name}[{i}]")

    if with_bounds:

        if not isneginf(var.lb):

            ax.plot(xValues,
                    [var.lb for _ in xValues],
                    color="orange",
                    label="lb")

        if not isinf(var.ub):

            ax.plot(xValues,
                    [var.ub for _ in xValues],
                    color="red",
                    label="ub")

    if with_legend:
        ax.legend()

    if with_labels:

        ax.set_xlabel(xlabel)

        if ylabel is None:

            ax.set_ylabel(var.name)

    if with_grid:
        ax.grid()

    if save:
        plt.savefig(fname=join(path, var.name))

    if show:
        plt.show()

    return ax
