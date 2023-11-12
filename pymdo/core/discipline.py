from typing import Dict
from typing import List
from typing import Tuple
from copy import deepcopy

from numpy import ndarray
from numpy import zeros
from numpy import imag

from .variable import Variable
from .variable import FLOAT_DATA_TYPE
from .variable import COMPLEX_FLOAT_DATA_TYPE
from .cache import cache_factory


class VariableMissingError(Exception):
    def __init__(self,
                 var: Variable,
                 discName: str,
                 isInput: bool = True
                 ) -> None:
        self.var: Variable = var
        self.varType: str = "input" if isInput else "output"
        self.message = f"In {discName}, {self.varType}: {self.var.name} is missing"
        super().__init__(self.message)


class VariableSizeError(Exception):
    def __init__(self,
                 var: Variable,
                 discName: str,
                 wrongSize: int,
                 isInput: bool = True
                 ) -> None:
        self.var: Variable = var
        self.varType: str = "input" if isInput else "output"
        self.message = f"In {discName}, {self.varType}: {self.var.name} (size: {self.var.size}) is passed with wrong size {wrongSize}"
        super().__init__(self.message)


class JacobianEntryError(Exception):
    def __init__(self,
                 outVar: Variable,
                 discName: str,
                 inVar: Variable = Variable("_"),
                 wrongSize: Tuple[int, int] = None
                 ) -> None:
        if wrongSize == None:
            self.message = f"In {discName}, jacobian entry: ({outVar}, {inVar}) is missing"
        else:
            self.message = f"In {discName}, jacobian entry: ({outVar}, {inVar}) (size: {outVar.size}, {inVar.size}) has wrong size {wrongSize}"
        super().__init__(self.message)


class Discipline:
    """

    Base discipline class

    """

    ANALYTIC = "Analytic"

    FINITE_DIFFERENCE = "FiniteDifference"

    COMPLEX_STEP = "ComplexStep"

    DIFF_METHODS = [ANALYTIC,
                    FINITE_DIFFERENCE,
                    COMPLEX_STEP]

    def __init__(self,
                 name: str,
                 inputVars: List[Variable],
                 outputVars: List[Variable],
                 cacheType: str = "memory",
                 cachePolicy: str = "latest",
                 cacheTol: float = 1e-9,
                 cachePath: str = None):

        self.name: str = name
        """ Discipline name by which it is accessed """

        self.nInputs: int = len(inputVars)
        """ Number of input variables """

        self.sizeInputs: int = sum([var.size for var in inputVars])
        """ Total size of input variables """

        self.inputVars: List[Variable] = inputVars
        """ List of input variables """

        self.nOutputs: int = len(outputVars)
        """ Number of output variables """

        self.sizeOutputs: int = sum([var.size for var in outputVars])
        """ Total size of input variables """

        self.outputVars: List[Variable] = outputVars
        """ List of discipline output variables """

        self.diffInputs: List[Variable] = []
        """ List of variables w.r.t differentiate """

        self.diffOutputs: List[Variable] = []
        """ List of variables to be differentiated """

        self.diffMethod: str = self.ANALYTIC
        """ Jacobian computation method """

        self.eps: float = 1e-4
        """ Jacobian approximation step (if needed) """

        self._floatDataType = FLOAT_DATA_TYPE
        """ Data type for floating point operations """

        self.defaultInputs: Dict[str, ndarray] = {}
        """ Default input values for evaluation """

        self.values: Dict[str, ndarray] = {}
        """ Latest evaluation input and output values """

        self.jac: Dict[str, Dict[str, ndarray]] = {}
        """ Latest evaluation jacobian """

        self.nEval: int = 0
        """ Number of calls to _eval() """

        self.nDiff: int = 0
        """ Number of calls to _differentiate() or _approximate_jacobian()"""

        if cachePath is None:
            cachePath = self.name

        self.cache = cache_factory(self.inputVars, 
                                   self.outputVars,
                                   cacheType,
                                   cachePolicy,
                                   cacheTol,
                                   cachePath)
        """ Discipline cache, stores evaluation values and jacobians """

        self.approximating: bool = False
        """ Flag that is turned on when _approximate_jacobian() is executing """

    def __repr__(self) -> str:
        return self.name

    def verify_values(self, vars: List[Variable]) -> None:
        """ Verify that values contains entries with correct size for all variables """

        for var in vars:

            if var.name not in self.values:
                isInput = var in self.inputVars
                raise VariableMissingError(var,
                                           self.name,
                                           isInput)

            if var.size != self.values[var.name].size:
                isInput = var in self.inputVars
                raise VariableSizeError(var,
                                        self.name,
                                        self.values[var.name].size,
                                        isInput)
            
    def verify_jac(self) -> None:
        """ Verify that the current jacobian contains entries with correct size for all diff. variables """

        for outVar in self.diffOutputs:

            if outVar.name not in self.jac:

                raise JacobianEntryError(outVar,
                                         self.name)

            for inVar in self.diffInputs:

                if inVar.name not in self.jac[outVar.name]:

                    raise JacobianEntryError(outVar,
                                             self.name,
                                             inVar)
                
                if self.jac[outVar.name][inVar.name].shape != (outVar.size, inVar.size):

                    raise JacobianEntryError(outVar,
                                             self.name,
                                             inVar,
                                             self.jac[outVar.name][inVar.name].shape)

    def load_cache(self, path: str = None) -> None:

        if self.cache is not None:

            self.cache.from_disk(path)
    
    def save_cache(self, path: str = None) -> None:

        if self.cache is not None:

            self.cache.to_disk(path)
            
    def get_input_values(self) -> Dict[str, ndarray]:
        if not self.values:
            return {}
        return {var.name: self.values[var.name] for var in self.inputVars}

    def get_output_values(self) -> Dict[str, ndarray]:
        if not self.values:
            return {}
        return {var.name: self.values[var.name] for var in self.outputVars}

    def get_input_names(self) -> List[str]:

        return [var.name for var in self.inputVars]

    def get_output_names(self) -> List[str]:

        return [var.name for var in self.outputVars]
    
    def get_diff_input_names(self) -> List[str]:

        return [var.name for var in self.diffInputs]

    def get_diff_output_names(self) -> List[str]:

        return [var.name for var in self.diffOutputs]
    
    def set_default_inputs(self, values: Dict[str, ndarray]) -> None:

        self.defaultInputs.update(deepcopy(values))

    def _load_inputs(self, inputValues: Dict[str, ndarray] = None):
        """
        Load the input values provided.

        If no value is provided for a variable,
        try to use the default value.

        Check that no values are missing and that the sizes
        are correct.

        """

        _inputValues = self.defaultInputs
        
        if inputValues is not None:

            _inputValues.update(inputValues)

        self.values.update(deepcopy(_inputValues))
        
        self.verify_values(self.inputVars)
        
    def _eval(self) -> None:
        raise NotImplementedError

    def eval(self, inputValues: Dict[str, ndarray] = None) -> Dict[str, ndarray]:
        """

        Execute discipline with the given inputValues.

        """

        self._load_inputs(inputValues)

        if self.cache is not None:

            outputValues, _ = self.cache.load_entry(self.get_input_values())

            if outputValues:
                
                self.values.update(outputValues)

                return self.values

        self._eval()

        self.verify_values(self.outputVars)
        
        self.nEval += 1

        if self.approximating == False:
            
            self.set_default_inputs(self.get_input_values())

            if self.cache is not None:
                
                self.cache.add_entry(self.get_input_values(),
                                 self.get_output_values(),
                                 None)
        return self.values

    def add_diff_inputs(self, vars: List[Variable] = None) -> None:
        """

        Add variable(s) w.r.t which to differentiate.

        If no variables are provided,
        all input variables are set as such


        """

        if vars is None:

            vars = self.inputVars
        
        inputNames = self.get_input_names()

        diffInputNames = self.get_diff_input_names()

        for var in vars:

            if var.name not in inputNames:

                raise VariableMissingError(var,
                                           self.name,
                                           True)
            
            if var.name not in diffInputNames:

                self.diffInputs.append(var)
        
        if self.cache is not None:

            self.cache.update_diff_inputs(self.diffInputs)

    def add_diff_outputs(self, vars: List[Variable] = None) -> None:
        """

        Add variable(s) to be differentiated.

        If no variables are provived,
        all output variables are set as such.

        """

        if vars is None:

            vars = self.outputVars
        
        outputNames = self.get_output_names()
        
        diffOutputNames = self.get_diff_output_names()

        for var in vars:

            if var.name not in outputNames:

                raise VariableMissingError(var,
                                           self.name,
                                           False)
            
            if var.name not in diffOutputNames:

                self.diffOutputs.append(var)
        
        if self.cache is not None:

            self.cache.update_diff_outputs(self.diffOutputs)

    def set_jacobian_approximation(self, method: str = FINITE_DIFFERENCE, eps: float = 1e-4) -> None:
        self.diffMethod = method
        self.eps = eps

    def _approximate_jacobian(self) -> None:

        F = {outVar.name: self.values[outVar.name]
             for outVar in self.outputVars}

        X = {inVar.name: self.values[inVar.name] for inVar in self.inputVars}

        if self.diffMethod == self.COMPLEX_STEP:

            self._floatDataType = COMPLEX_FLOAT_DATA_TYPE

            X = {name: values.astype(self._floatDataType)
                 for name, values in X.items()}

            Eps = self.eps * 1j

        else:

            Eps = self.eps

        for inVar in self.diffInputs:

            for i in range(0, inVar.size):


                Xp = deepcopy(X)

                Xp[inVar.name][i] += Eps

                Fp = self.eval(Xp)

                for outVar in self.diffOutputs:

                    if self.diffMethod == self.COMPLEX_STEP:

                        dF = imag(Fp[outVar.name])

                    else:

                        dF = Fp[outVar.name] - F[outVar.name]

                    self.jac[outVar.name][inVar.name][:, i] = dF / self.eps

        if self.diffMethod == self.COMPLEX_STEP:

            self._floatDataType = FLOAT_DATA_TYPE

            X = {name: values.real
                 for name, values in X.items()}

        # Restore self.values
        self.values.update(X)
        self.values.update(F)

    def _init_jacobian(self) -> None:
        """

        If differentiated inputs/outputs are not defined,
        all input/output variables are set as such.

        """

        if not self.diffInputs:
            self.add_diff_inputs()

        if not self.diffOutputs:
            self.add_diff_outputs()

        self.jac = {outVar.name:
                    {inVar.name:
                     zeros(
                         (outVar.size, inVar.size), self._floatDataType)
                     for inVar in self.diffInputs}
                    for outVar in self.diffOutputs}

        return self.jac

    def _differentiate(self) -> None:
        raise NotImplementedError

    def differentiate(self, inputValues: Dict[str, ndarray] = None) -> Dict[str, Dict[str, ndarray]]:
        """

        Differentiate the discipline for a given set of input values.

        If no input values are provided, try to use the inputs
        from the latest evaluation.

        The discipline must be evaluated for the given values.

        """

        # Avoid re-evaluation at latest point
        if not (inputValues is None and self.values):

            self.eval(inputValues)

        if self.cache is not None:

            _, jac = self.cache.load_entry(self.get_input_values())

            if jac:

                self.jac = jac

                return self.jac

        self.jac = self._init_jacobian()

        if self.diffMethod == self.ANALYTIC:
            self._differentiate()
        else:
            self.approximating = True
            self._approximate_jacobian()
            self.approximating = False
            
        self.verify_jac()

        self.nDiff += 1

        if self.cache is not None:

            self.cache.add_entry(self.get_input_values(),
                                 None,
                                 self.jac)

        return self.jac