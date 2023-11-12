from typing import Dict
from typing import List

from numpy import ndarray

from pymdo.core.variable import Variable
from pymdo.core.discipline import Discipline
from pymdo.mda.mda import MDA
from pymdo.mda.smart_mda import SmartMDA
from .mdo import MDOptProblem


class MDF(MDOptProblem):

    def __init__(self,
                 disciplines: List[Discipline],
                 designVariables: List[Variable],
                 objective: Variable,
                 maximizeObjective: bool = False,
                 useNormalization: bool = True,
                 saveDesignVector: bool = False) -> None:

        super().__init__(disciplines,
                         designVariables,
                         objective,
                         maximizeObjective,
                         useNormalization,
                         saveDesignVector)

        self.mda: MDA = SmartMDA(self.disciplines)
        """ MDA used for converging the system to feasibility """

    def _set_values(self) -> Dict[str, ndarray]:

        return self.mda.eval(self.designVector)
    
    def _set_grad(self) -> Dict[str, Dict[str, ndarray]]:

        self.mda.add_diff_inputs(self.designVariables)

        self.mda.add_diff_outputs([self.objective])

        self.mda.add_diff_outputs(self.constraints)
        
        return self.mda.differentiate()