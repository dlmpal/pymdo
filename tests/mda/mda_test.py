import numpy as np 
from pymdo.core.variable import Variable
from pymdo.core.discipline import Discipline
from pymdo.mda.factory import mda_factory
from pymdo.mda.smart_mda import SmartMDA
import matplotlib.pyplot as plt

plt.style.use("dark_background")

def func1(input_vars):
    x1 = input_vars["x1"][0]
    z = input_vars["z"][0]
    y12 = input_vars["y12"][0]
    y21 = x1**2 + x1*z - y12 * z
    return {"y21": y21}


def dfunc1(input_vars):
    x1 = input_vars["x1"][0]
    z = input_vars["z"][0]
    y12 = input_vars["y12"][0]

    dx1 = 2 * x1 + z
    dz = np.array([x1 - y12])
    dy12 = np.array([-z])
    jac = {}
    jac["y21"] = {"x1": dx1, "z": dz, "y12": dy12}

    return jac


def func2(input_vars):

    x2 = input_vars["x2"][0]
    z = input_vars["z"][0]
    y21 = input_vars["y21"][0]
    y12 = 2*y21 - x2**2 + z*x2
    return {"y12": np.array([y12])}


def dfunc2(input_vars):
    x2 = input_vars["x2"][0]
    z = input_vars["z"][0]
    y21 = input_vars["y21"][0]
    dx2 = np.array([-2*x2 + z])
    dz = np.array([x2])
    dy21 = np.array([2])
    jac = {}
    jac["y12"] = {"x2": dx2, "z": dz, "y21": dy21}
    return jac

class Disc1(Discipline):

    def __init__(self, _inputVars, _outputVars):

        super().__init__("D1",
                         _inputVars,
                         _outputVars,
                         cacheType="memory",
                         cachePolicy="full")

    def _eval(self) -> None:
        self.values.update(func1(self.values))

    def _differentiate(self) -> None:
        self.jac.update(dfunc1(self.values))


class Disc2(Discipline):

    def __init__(self, _inputVars, _outputVars):

        super().__init__("D2",
                         _inputVars,
                         _outputVars,
                         cacheType="memory",
                         cachePolicy="full")

    def _eval(self) -> None:
        self.values.update(func2(self.values))

    def _differentiate(self) -> None:
        self.jac.update(dfunc2(self.values))


def ConfirmTotalDerivatives(input_vars):
    x1 = input_vars["x1"][0]
    z = input_vars["z"][0]
    x2 = input_vars["x2"][0]
    y12 = input_vars["y12"][0]
    y21 = input_vars["y21"][0]

    # Total
    dY12dx1 = (4*x1+2*z)/(1 + 2*z)
    dY12dx2 = (-2*x2+z)/(1+2*z)
    dY21dx1 = (2*x1+z)/(2*z+1)
    dY21dx2 = (-z**2 + 2*z*x2)/(2*z+1)
    dY12dz = (2*x1+x2-4*x1**2+2*x2**2)/(1+2*z)**2
    dY21dz = (x2**2 - 2*x2*z+x1-2*x1**2 - 2*z**2 * x2)/(1+2*z)**2

    # dRdY
    dy21dy12 = -z
    dy12dy21 = 2.0
    dRdY = np.array([[1.0, -dy21dy12], [-dy12dy21, 1.0]], np.float64)
    print(dRdY)

    # dRdX
    dx1 = 2*x1 + z
    dz1 = x1 - y12
    dx2 = -2*x2 + z
    dz2 = x2
    dRdX = np.array([[dx1, 0.0, dz1],
                     [0.0, dx2, dz2]])
    print(dRdX)
    
    # dFdY
    dFdY = np.array([[1.0, 0.0],
                     [0.0, 1.0]])
    print(dFdY)
    dFdX = dFdY @ np.linalg.solve(dRdY, dRdX)

    print(dFdX)

    return {"y12": {"x1": dY12dx1, "x2": dY12dx2, "z": dY12dz}}, {"y21": {"x1": dY21dx1, "x2": dY21dx2, "z": dY21dz}}


x1 = Variable("x1", 1)
x2 = Variable("x2", 1)
z = Variable("z", 1)
y12 = Variable("y12", 1)
y21 = Variable("y21", 1)

disc1 = Disc1([x1, z, y12], [y21])
disc2 = Disc2([x2, z, y21], [y12])

inputs = {"x1": np.array([3.0]), "x2": np.array([3.0]), "z": np.array(
    [3.0]), "y12": np.array([1.0]), "y21": np.array([1.0])}

disc1.set_jacobian_approximation(Discipline.FINITE_DIFFERENCE)

disc2.set_jacobian_approximation(Discipline.COMPLEX_STEP)

mda = mda_factory([disc2, disc1], mdaType="MDAGaussSeidel", saveVariableLog = True)

mda.set_options(150,
               0.1,
               1e-6)

mda.eval(inputs)

mda.plot_residual(["total"], False, True)

mda.plot_variable(False, True)

disc1.save_cache()
disc2.save_cache()

print(disc1.nEval, disc1.nDiff, disc2.nEval, disc2.nDiff)