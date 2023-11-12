from typing import List, Dict 

from numpy import ndarray

from pymdo.core.discipline import Discipline
from .mda import MDA


class MDAGaussSeidel(MDA):
    """
    
    This MDA sub-class implements the generalized 
    or non-linear Gauss-Seidel iteration:

    Yi^(k+1) = Yi(Xi^k),

    where Xi^k = [xi^(k+1) z^(k+1) y1i^(k+1) ... y(i-1)i^(k+1)  y(i+1)i^k yni^k]

    """

    def __init__(self,
                 disciplines: List[Discipline],
                 name: str = "MDAGaussSeidel",
                 nIterMax: int = 15,
                 relaxFact: float = 1.0,
                 relTol: float = 0.0001,
                 saveVariableLog: bool = False) -> None:

        super().__init__(disciplines,
                         name,
                         nIterMax,
                         relaxFact,
                         relTol,
                         saveVariableLog)

    def _eval(self) -> None:

        self.mdaStatus = self.MDA_FAIL

        self.residualLog = []

        self.variableLog = []
        
        while self._terminate_condition() == False:

            currentOutputs: Dict[str, ndarray] = {}

            for disc in self.disciplines:

                discInputs = {var.name: self.values[var.name] if var.name not in currentOutputs
                              else currentOutputs[var.name] for var in disc.inputVars}

                discOutputs = disc.eval(discInputs)

                for var in disc.outputVars:

                    currentOutputs[var.name] = self.relaxFact * discOutputs[var.name] + \
                        (1 - self.relaxFact) * self.values[var.name]

            self._compute_residual(currentOutputs,
                                  self.values)

            self.values.update(currentOutputs)