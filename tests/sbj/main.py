from typing import Dict, List, Callable

import numpy as np
import matplotlib.pyplot as plt
from pymdo.core.discipline import Discipline
from pymdo.core.cache import cache_factory
from pymdo.mda.smart_mda import SmartMDA
from pymdo.mda.factory import mda_factory
from pymdo.mdo.mdf import MDF

from functions import *


plt.style.use("dark_background")

class SBJ_Discipline(Discipline):

    def __init__(self,
                 _name: str,
                 _inputNames: List[str],
                 _outputNames: List[str],
                 _func: Callable[[Dict[str, np.ndarray]], Dict[str, np.ndarray]]):

        inputVars = [Variable(name, 1) for name in _inputNames]
        outputVars = [Variable(name, 1) for name in _outputNames]

        self.func: Callable[[Dict[str, np.ndarray]],
                            Dict[str, np.ndarray]] = _func

        super().__init__(_name,
                         inputVars,
                         outputVars)

    def _eval(self) -> None:
        self.values.update(self.func(self.values))


flightConditions = SBJ_Discipline("flightConditions",
                                  [z_cr_str, M_cr_str],
                                  [T_cr_str, p_cr_str, rho_cr_str, V_cr_str],
                                  ExecuteFlightConditions)

wingGeometry = SBJ_Discipline("wingGeometry",
                              [S_w_str, fwLE_w_str, fwTE_w_str, tr_w_str,
                               tcr_w_str],
                              [c_w_str, b_w_str, AR_w_str,
                                  fw25_w_str, Sexp_w_str, Swet_w_str],
                              ExecuteWingGeometry)

vertStabGeometry = SBJ_Discipline("vertStabGeometry",
                                  [S_w_str, fwLE_v_str, fwTE_v_str, tcr_v_str,
                                   tr_v_str],
                                  [S_v_str, c_v_str, b_v_str,
                                      Swet_v_str, fw25_v_str],
                                  ExecuteVerticalStabilizerGeometry)


fuselageGeometry = SBJ_Discipline("fuselageGeometry",
                                  [W_fuel_str, M_cr_str],
                                  [L_f_str, S_f_str],
                                  ExecuteFuselageGeometry
                                  )

weight = SBJ_Discipline("weight",
                        [S_w_str, tcr_w_str, fw25_w_str, Sexp_w_str, tr_w_str, b_w_str,
                         S_v_str, b_v_str, tcr_v_str, fw25_v_str,
                         S_f_str, L_f_str,
                         p_cr_str,
                         P0_str,
                         W_fuel_str],
                        [tow_str, zfw_str],
                        ExecuteWeight
                        )

weight.set_default_inputs({P0_str: np.array([10e3])})

engine = SBJ_Discipline("engine",
                        [p_cr_str, rho_cr_str, T_cr_str,
                         M_cr_str, S_w_str, V_cr_str, CD_str],
                        [sfc_str, P_str, P0_str, L_nac_str, D_nac_str, fa_eng_str],
                        ExecuteEngine
                        )

engine.set_default_inputs({CD_str: np.array([0.02]),
                             sfc_str: np.array([0.01]),
                             P_str: np.array([9e3])})

aerodynamics = SBJ_Discipline("aerodynamics", [S_w_str, Swet_w_str, c_w_str, tcr_w_str, b_w_str, tr_w_str, fwLE_w_str, fwTE_w_str, AR_w_str,
                                              S_v_str, Swet_v_str, c_v_str, tcr_v_str, b_v_str, tr_v_str, fwLE_v_str, fwTE_v_str,
                                              L_f_str, S_f_str,
                                              L_nac_str, D_nac_str, fa_eng_str,
                                              M_cr_str, T_cr_str, rho_cr_str, V_cr_str,
                                              tow_str],
                              [CL_str, CD_str],
                              ExecuteAerodynaics
                              )

aerodynamics.set_default_inputs({L_nac_str: np.array([3.0]),
                                   D_nac_str: np.array([1.0]),
                                   fa_eng_str: np.array([1.0]),
                                   tow_str: np.array([30e3])})

performance = SBJ_Discipline("performance",
                             [V_cr_str, sfc_str, CL_str,
                                 CD_str, tow_str, W_fuel_str],
                             [Range_str],
                             ExecutePerformance)

disciplines = [flightConditions, wingGeometry, vertStabGeometry,
               fuselageGeometry, engine, weight, aerodynamics, performance]

for disc in disciplines:

    disc.set_jacobian_approximation("ComplexStep", 1e-3)

mode = 1 # 0 for mda, 1 for mdo #

if mode == 0:

    initialDesign, _ = create_design_vars()

    mdaList = {
        "MDAJacobi": {"relaxFact": 1.0, "nIterMax": 10, "relTol": 1e-3 * 0.5},
        "MDANewton": {"nIterMax": 10, "relTol": 1e-3},
        "MDAGaussSeidel": {"relaxFact": 1.0, "nIterMax": 10, "relTol": 1e-3 * 0.5}}

    mdaLogs = {}

    for name, kwargs in mdaList.items():

        mda = SmartMDA(disciplines)

        mda.groups[1]["Group_1"] = mda_factory([engine, aerodynamics, weight], name, **kwargs)

        inputs = {}

        inputs.update(initialDesign)

        mda.eval(inputs)

        print(aerodynamics.nEval)

        mdaLogs.update({name: mda.groups[1]["Group_1"].residualLog})

    for name, log in mdaLogs.items():

        iters = [i for i in range(len(log))]

        total = [log[i]["CD"] for i in iters]

        plt.semilogy(iters, total, label=name)

    plt.xlabel("iterations")
    plt.ylabel("Drag coeff residual")
    plt.legend()
    plt.show()

if mode == 1:

    Range = Variable(Range_str, 1)

    TOW = Variable(tow_str, 1, ub = 49e3)

    initialDesign, designVars = create_design_vars()

    mdf = MDF(disciplines,
              designVars,
              Range,
              True,
              True,
              True)
    
    mdf.mda = mda_factory(mdf.disciplines, "MDAGaussSeidel", relaxFact = 0.9)
    
    mdf.mda.cache = cache_factory(mdf.mda.inputVars,
                                  mdf.mda.outputVars,
                                  "memory",
                                  "full",
                                  path = "cache/mda")
    mdf.add_constraint(TOW)

    mdf.execute(initialDesign, "slsqp", options={"maxiter": 5})
    
    mdf.mda.save_cache()
    
    mdf.plot_optimization_history(show = False, save=True, path = "results")
