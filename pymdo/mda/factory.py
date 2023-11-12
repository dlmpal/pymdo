from typing import List

from pymdo.core.discipline import Discipline
from .mda import MDA 
from .gauss_seidel import MDAGaussSeidel
from .jacobi import MDAJacobi
from .newton import MDANewton
from .hybrid import MDAHybrid

class InvalidMDAName(Exception):

    def __init__(self, _invalidMDAType: str) -> None:

        self.message = f"Invalid MDA type: {_invalidMDAType}. Available types are: MDAGaussSeidel, MDAJacobi, MDANewton, MDAHybrid"

        super().__init__(self.message)
        
def mda_factory(disciplines: List[Discipline], mdaType: str = "MDAGaussSeidel", **kwargs) -> MDA:
    
    if mdaType == "MDAGaussSeidel":

        return MDAGaussSeidel(disciplines, **kwargs)
    
    if mdaType == "MDAJacobi":

        return MDAJacobi(disciplines, **kwargs)
    
    if mdaType == "MDANewton":

        return MDANewton(disciplines, **kwargs)
    
    if mdaType == "MDAHybrid":

        return MDAHybrid(disciplines, **kwargs)

    else:
        raise InvalidMDAName(mdaType)