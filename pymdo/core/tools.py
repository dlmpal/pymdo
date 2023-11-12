from typing import Dict, List, Tuple

from .variable import Variable
from .discipline import Discipline


def seperate_inputs_and_outputs(_disciplines: List[Discipline]) -> Tuple[List[Variable], List[Variable]]:
    """ 

    Find input and output variables from list of disciplines

    """
    outputVarDict: Dict[str, Variable] = {}
    for disc in _disciplines:
        for var in disc.outputVars:
            outputVarDict[var.name] = var

    inputVarDict: Dict[str, Variable] = {}
    for disc in _disciplines:
        for var in disc.inputVars:
            if var.name not in outputVarDict:
                inputVarDict[var.name] = var

    return (inputVarDict.values(), outputVarDict.values())