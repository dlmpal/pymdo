from typing import List
from typing import Dict

from numpy import ndarray
from numpy import zeros
from numpy import reshape
from scipy.sparse import dok_matrix # To be changed to dok_array
from scipy.sparse import eye
from scipy.sparse.linalg import spsolve

from .variable import Variable 
from .variable import FLOAT_DATA_TYPE
from .variable import array2d_to_dict
from .discipline import Discipline


class DerivativeAssembler():
    """

    Assemble derivative matrices for a set of disciplines.

    This object has two main usecases:
    
    1) Computing the coupled derivatives for MDA objects (dFdX)

    2) Assembling the Jacobian matrix for MDANewton

    If the same instance is to be used multiple times, use update_jac
    before calling the rest of the methods, so that the discipline jacobians 
    are up to date.

    If the differentiated inputs/outputs or coupling variables change, 
    a new instance must be created, or the results will be incorrect.

    """

    def __init__(self,
                 disciplines: List[Discipline],
                 couplings: List[Variable],
                 diffInputs: List[Variable],
                 diffOutputs: List[Variable]):

        self.disciplines = disciplines
        """ List of disciplines for which the coupled derivatives are to computed """

        self.couplings = couplings
        """ Coupling variables """

        self.Ny = sum([y.size for y in self.couplings])
        """ Total size of all couplings """

        self.diffInputs = diffInputs
        """ Variables w.r.t differentiate """

        self.Nx = sum([x.size for x in self.diffInputs])
        """ Total size of differentiated inputs """

        self.diffOutputs = diffOutputs
        """ Variables to be differentiated """

        self.Nf = sum([f.size for f in self.diffOutputs])
        """ Total size of differentiated outputs """

        self.jac: Dict[str, Dict[str, ndarray]] = self.update_jac()
        """ Partial derivatives of the couplings 
            w.r.t other couplings and diffInputs
        """

    def update_jac(self) -> Dict[str, Dict[str, ndarray]]:

        self.jac = {}

        for disc in self.disciplines:

            self.jac.update(disc.jac)
        
        return self.jac
    
    def dYdY(self) -> dok_matrix:
        """ 

        dR/dY

        """

        _dYdY = dok_matrix((self.Ny,
                           self.Ny),
                          dtype=FLOAT_DATA_TYPE)

        r = 0

        for yi in self.couplings:

            c = 0

            for yj in self.couplings:

                if yi == yj:

                    _dYdY[r: r + yi.size, c: c +
                          yj.size] = eye(yi.size, dtype=FLOAT_DATA_TYPE)

                else:

                    if yj.name in self.jac[yi.name]:

                        _dYdY[r: r + yi.size, c: c +
                              yj.size] = -self.jac[yi.name][yj.name]

                c += yj.size

            r += yi.size

        return _dYdY

    def dYdX(self) -> dok_matrix:
        """ 

        dR/dX 

        """

        _dYdX = dok_matrix((self.Ny,
                           self.Nx),
                          dtype=FLOAT_DATA_TYPE)

        r = 0

        for yi in self.couplings:

            c = 0

            for xj in self.diffInputs:
                
                if xj.name in self.jac[yi.name]:

                    _dYdX[r: r + yi.size, c: c +
                          xj.size] = self.jac[yi.name][xj.name]

                c += xj.size
            
            r += yi.size

        return _dYdX

    def dFdY(self) -> dok_matrix:
        """ 

        dF/dY

        """

        _dFdY = dok_matrix((self.Nf,
                           self.Ny),
                          dtype=FLOAT_DATA_TYPE)
        r = 0

        for fi in self.diffOutputs:

            c = 0

            for yj in self.couplings:

                if fi == yj:

                    _dFdY[r: r + fi.size, c: c + yj.size] = eye(fi.size, dtype=FLOAT_DATA_TYPE)

                c += yj.size
        
            r += fi.size

        return _dFdY

    def dFdX(self):
        """ 

        dF/dX 

        """

        _dFdX = zeros((self.Nf,
                       self.Nx),
                      FLOAT_DATA_TYPE)

        _dYdY = self.dYdY()

        _dYdX = self.dYdX()

        _dFdY = self.dFdY().tocsr()

        if self.Nx >= self.Nf:
            """ Adjoint approach """
            _dYdY = _dYdY.tocsr()
            _dYdX = _dYdX.tocsr()
            _dFdX = spsolve(_dYdY.transpose(), _dFdY.transpose()).transpose() @ _dYdX

        else:
            """ Direct approach"""
            _dYdY = _dYdY.tocsc()
            _dYdX = _dYdX.tocsc()
            _dFdX = _dFdY @ spsolve(_dYdY, _dYdX)
        
        if type(_dFdX) != ndarray:

            _dFdX = _dFdX.toarray()

        _dFdX = reshape(_dFdX, (self.Nf, self.Nx))
        
        dFdX_Dict = array2d_to_dict(self.diffInputs,
                                 self.diffOutputs,
                                 _dFdX)

        return dFdX_Dict
