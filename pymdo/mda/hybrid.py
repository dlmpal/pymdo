from typing import List

from pymdo.core.discipline import Discipline
from .mda import MDA
from .gauss_seidel import MDAGaussSeidel
from .newton import MDANewton


class MDAHybrid(MDA):
    """

    A hybrid MDA consisting of a sequence of MDAs.

    By default a GaussSeidel MDA is followed by a
    Newton MDA. This combination uses the robust nature 
    of the former, with the fast convergence of the latter.

    Any user-defined sequence of MDAs can be used.

    """

    def __init__(self,
                 disciplines: List[Discipline],
                 name: str = "MDAHybrid",
                 mdaSequence: List[MDA] = None,
                 saveVariableLog: bool = False) -> None:

        super().__init__(disciplines,
                         name,
                         saveVariableLog=saveVariableLog)

        self.mdaSequence: List[MDA] = mdaSequence

        if mdaSequence is None:

            self.mdaSequence = [MDAGaussSeidel(self.disciplines,
                                               nIterMax=2,
                                               relaxFact=0.8,
                                               relTol=10,
                                               saveVariableLog=saveVariableLog),

                                MDANewton(self.disciplines,
                                          saveVariableLog=saveVariableLog)]

    def _eval(self) -> None:

        self.residualLog = []

        self.variableLog = []

        for mda in self.mdaSequence:

            self.values.update(mda.eval(self.values))

            self.residualLog.extend(mda.residualLog)

            if self.saveVariableLog:

                self.variableLog.extend(mda.variableLog)
