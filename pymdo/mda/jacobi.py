from typing import List 
from typing import Dict 

from numpy import ndarray

from pymdo.core.discipline import Discipline
from .mda import MDA

class MDAJacobi(MDA):
    """
    
    This MDA sub-class implements the generalized 
    or non-linear Jacobi iteration:

    Yi^(k+1) = Yi(Xi^k),

    where Xi^k = [xi^k z^k y1i^k ... yni^k], j =/= i

    """

    def __init__(self,
                 disciplines: List[Discipline],
                 name: str = "MDAJacobi",
                 nIterMax: int = 15,
                 relaxFact: float = 1.0,
                 relTol:float = 0.0001,
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

                discInputs = {
                    var.name: self.values[var.name] for var in disc.inputVars}

                disc.eval(discInputs)
                
                currentOutputs.update({var.name: disc.values[var.name] for var in disc.outputVars})
                
            for var in self.outputVars:
                currentOutputs[var.name] = self.relaxFact * currentOutputs[var.name] + \
                    (1 - self.relaxFact) * self.values[var.name]
            
            self._compute_residual(currentOutputs,
                                  self.values)

            self.values.update(currentOutputs)
