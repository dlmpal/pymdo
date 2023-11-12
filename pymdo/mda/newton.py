from typing import List
from typing import Dict

from numpy import ndarray
from scipy.sparse.linalg import spsolve

from pymdo.core.discipline import Discipline
from pymdo.core.derivative_assembler import DerivativeAssembler
from .mda import MDA


class MDANewton(MDA):
    """

    This MDA sub-class uses the Newton iteration 
    for a system of non-linear equations:

    dR^k/dY * Ycorr^(k) = R^k

    Y^(k+1) = Y^k + Ycorr^k

    """

    def __init__(self,
                 disciplines: List[Discipline],
                 name: str = "MDANewton",
                 nIterMax: int = 15,
                 relTol: float = 0.0001,
                 saveVariableLog: bool = False) -> None:

        super().__init__(disciplines,
                         name,
                         nIterMax,
                         relTol=relTol,
                         saveVariableLog=saveVariableLog)

    def _eval(self) -> None:

        self.mdaStatus = self.MDA_FAIL

        self.residualLog = []

        self.variableLog = []

        assembler = DerivativeAssembler(self.disciplines,
                                        self.outputVars,
                                        self.diffInputs,
                                        self.diffOutputs)

        while self._terminate_condition() == False:

            currentOutputs: Dict[str, ndarray] = {}

            for disc in self.disciplines:

                disc.eval(self.values)

                disc.differentiate()

                currentOutputs.update(disc.get_output_values())

            R = self._compute_residual(currentOutputs,
                                       self.values)

            assembler.update_jac()

            dRdY = assembler.dYdY().tocsr()

            Ycorr = spsolve(dRdY, R)

            r = 0

            for var in self.outputVars:

                self.values[var.name] = self.values[var.name] + \
                    Ycorr[r: r + var.size]

                r += var.size
