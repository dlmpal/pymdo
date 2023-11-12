from typing import Dict
from typing import List
from typing import Tuple
from typing import Callable

from numpy import ndarray
from numpy import zeros
from numpy import ones
from numpy import bool8
from scipy.optimize import Bounds
from scipy.optimize import NonlinearConstraint
from scipy.optimize import minimize

from pymdo.core.variable import Variable
from pymdo.core.variable import FLOAT_DATA_TYPE
from pymdo.core.variable import array_to_dict
from pymdo.core.variable import dict_to_array
from pymdo.core.variable import dict_to_array2d
from pymdo.core.variable import normalize_design_vector
from pymdo.core.variable import denormalize_design_vector
from pymdo.core.variable import normalize_gradient
from pymdo.core.variable import single_variable_plot
from pymdo.core.discipline import Discipline


class MDOptProblem:

    def __init__(self,
                 disciplines: List[Discipline],
                 designVariables: List[Variable],
                 objective: Variable,
                 maximizeObjective: bool = False,
                 useNormalization: bool = True,
                 saveDesignVector: bool = False) -> None:

        self.disciplines: List[Discipline] = disciplines
        """ List of disciplines modelling the problem """

        self.designVariables: List(Variable) = designVariables
        """ List of design variables """

        self.sizeDesignVars: int = sum(
            [dv.size for dv in self.designVariables])

        self.objective: Variable = objective
        """ Optimization objective """

        self.maximizeObjective: bool = maximizeObjective
        """ Whether to maximize the objective """

        self.constraints: List[Variable] = []
        """ List of constraints """

        self.useNormalization: bool = useNormalization
        """ Whether or not to normalize design variable and gradient values """

        self.values: Dict[str, ndarray] = {}
        """ Current values of discipline variables """

        self.grad: Dict[str, Dict[str, ndarray]] = {}
        """ Current gradient values """

        self.designVector: Dict[str, ndarray] = {}
        """ Current design variables values """

        self.optLog: List[Dict[str, ndarray]] = []
        """ Optimization Log.

            By default saves the objective and constraint values
            for each optimization cycle.

            Set saveDesignVector to True to also save the design vector.
        """

        self.saveDesignVector: bool = saveDesignVector
        """ Whether to save the design vector for each optimization cycle """

    def _get_design_variable_bounds(self) -> Bounds:

        lb = zeros(self.sizeDesignVars, FLOAT_DATA_TYPE)

        ub = zeros(self.sizeDesignVars, FLOAT_DATA_TYPE)

        keepFeasible = zeros(self.sizeDesignVars, bool8)

        r = 0

        for var in self.designVariables:

            _lb = 0.0 if self.useNormalization else var.lb
            _ub = 1.0 if self.useNormalization else var.ub

            lb[r: r + var.size] = _lb * ones(var.size, FLOAT_DATA_TYPE)

            ub[r: r + var.size] = _ub * ones(var.size, FLOAT_DATA_TYPE)

            keepFeasible[r: r + var.size] = var.keepFeasible * \
                ones(var.size, bool8)

            r += var.size

        return Bounds(lb,
                      ub,
                      keepFeasible)

    def add_constraint(self, constraint: Variable) -> None:
        """

        Add a constraint to the optimization problem.

        A constraint variable's behaviour is represented by its lower and upper bounds.

        If lb == ub, then it is treated as equality constraint (h(x) = 0).

        If lb = -inf and ub is finite, then the constraint is h(x) <= ub.

        Conversely, if ub = inf and lb is finite, then the constraint is h(x) >= lb

        """
        if constraint not in self.constraints:

            self.constraints.append(constraint)

    def _create_constraint_func(self, constraint: Variable) -> Tuple[Callable, Callable]:

        def h(inputValues: ndarray) -> ndarray:

            if self.values:

                return self.values[constraint.name]

            else:

                return zeros(constraint.size,
                             dtype=FLOAT_DATA_TYPE)

        def dh(inputValues: ndarray) -> ndarray:

            if self.grad:

                return dict_to_array2d(self.designVariables,
                                       [constraint],
                                       self.grad,
                                       True)

        return h, dh

    def _get_constraints(self) -> List[NonlinearConstraint]:

        cons = []

        for con in self.constraints:

            h, dh = self._create_constraint_func(con)

            _lb = con.lb

            _ub = con.ub

            cons.append(NonlinearConstraint(fun=h,
                                            lb=ones(con.size,
                                                    dtype=FLOAT_DATA_TYPE) * _lb,
                                            ub=ones(con.size,
                                                    dtype=FLOAT_DATA_TYPE) * _ub,
                                            jac=dh,
                                            keep_feasible=ones(con.size, bool8) * con.keepFeasible))

        return cons

    def _set_values(self) -> Dict[str, ndarray]:
        raise NotImplementedError

    def _set_grad(self) -> Dict[str, Dict[str, ndarray]]:
        raise NotImplementedError

    def _get_func(self) -> Tuple[Callable, Callable]:

        def F(designVectorArray: ndarray) -> float:

            self.designVector = array_to_dict(self.designVariables,
                                              designVectorArray)

            if self.useNormalization:
                self.designVector = denormalize_design_vector(self.designVariables,
                                                              self.designVector)

            self.values = self._set_values()

            return self._update_optimization_log()

        def dF(designVectorArray: ndarray) -> ndarray:

            self.grad = self._set_grad()

            if self.useNormalization:

                self.grad = normalize_gradient(self.designVariables,
                                               [self.objective] +
                                               self.constraints,
                                               self.grad)

            grad = dict_to_array2d(self.designVariables,
                                   [self.objective],
                                   self.grad,
                                   True)

            if self.maximizeObjective:

                grad = - grad

            return grad

        return F, dF

    def execute(self, initialDesignVector: Dict[str, ndarray], algoName: str = "SLSQP", options=None) -> Tuple[Dict[str, ndarray], float]:

        if self.useNormalization:
            initialDesignVector = normalize_design_vector(self.designVariables,
                                                          initialDesignVector)

        initialDesignVectorArray = dict_to_array(self.designVariables,
                                                 initialDesignVector)

        bnds = self._get_design_variable_bounds()

        cons = self._get_constraints()

        F, dF = self._get_func()

        result = minimize(fun=F,
                          x0=initialDesignVectorArray,
                          method=algoName,
                          jac=dF,
                          bounds=bnds,
                          constraints=cons,
                          options=options)

        return (self.designVector, result.fun)

    def _update_optimization_log(self) -> float:
        """
        Update the optimization log.

        Return the (signed) objecgive value.

        """

        self.optLog.append(
            {self.objective.name: self.values[self.objective.name]})

        self.optLog[-1].update({con.name: self.values[con.name]
                               for con in self.constraints})

        if self.saveDesignVector:
            self.optLog[-1].update(self.designVector)

        objValue = self.optLog[-1][self.objective.name]

        if self.maximizeObjective:

            objValue = - objValue

        return objValue

    def plot_optimization_history(self,
                                  show: bool = True,
                                  save: bool = False,
                                  path: str = "",
                                  with_bounds=True) -> None:

        for var in self.constraints + [self.objective]:

            values = {var.name: [self.optLog[i][var.name] for i in range(len(self.optLog))]
                      for var in self.constraints + [self.objective]}

            single_variable_plot(var,
                                 values[var.name],
                                 with_bounds=with_bounds,
                                 xlabel="cycles",
                                 show=show,
                                 save=save,
                                 path=path
                                 )

        if self.saveDesignVector:

            values = {var.name: [self.optLog[i][var.name] for i in range(len(self.optLog))]
                      for var in self.designVariables}

            for var in self.designVariables:

                single_variable_plot(var,
                                     values[var.name],
                                     with_bounds=with_bounds,
                                     xlabel="cycles",
                                     show=show,
                                     save=save,
                                     path=path
                                     )
