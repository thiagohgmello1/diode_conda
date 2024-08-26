import time

import numpy as np

from datetime import date
from model.optimizer import Optimizer
from model.material import Material
from model.particle import Particle

from pymoo.optimize import minimize
from pymoo.operators.mutation.pm import PM
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.termination import get_termination
from pymoo.visualization.scatter import Scatter
from pymoo.problems.functional import FunctionalProblem
from pymoo.operators.repair.rounding import RoundingRepair
from pymoo.operators.sampling.rnd import IntegerRandomSampling
from pymoo.operators.sampling.rnd import FloatRandomSampling


class MultiObjOpt(Optimizer):
    def __init__(self, params: dict, material: Material, particle_model: Particle, convergence: dict, scale=float):
        self.consts = self.constraints(params.pop("constraints"))
        self.boundaries = self.build_boundaries(params.pop("bounds"))
        super().__init__(
            **params, material=material, particle_model=particle_model, convergence=convergence, scale=scale
        )
        self.n_var = len(self.boundaries[0])
        self.obj_funcs = self.choose_objective_funcs()


    def choose_objective_funcs(self):
        obj_funcs = list()
        if "ZBI" in self.objectives["methods"]:
            obj_funcs.append(self.zero_voltage_imp)
        if "ZBR" in self.objectives["methods"]:
            obj_funcs.append(self.zero_bias_responsivity)
        if "ASY" in self.objectives["methods"]:
            obj_funcs.append(self.asymmetry)
        return obj_funcs


    def optimize(self):
        problem = FunctionalProblem(
            self.n_var,
            self.obj_funcs,
            constr_ieq=self.consts,
            xl=self.boundaries[0],
            xu=self.boundaries[1],
            callback=self.save_current_iter
        )
        algorithm = NSGA2(
            pop_size=self.pop_size,
            n_offsprings=10,
            sampling=IntegerRandomSampling(),
            crossover=SBX(prob=0.9, eta=15, vtype=float, repair=RoundingRepair()),
            mutation=PM(eta=20, vtype=float, repair=RoundingRepair()),
            eliminate_duplicates=True
        )
        termination = get_termination("n_gen", self.max_iter)

        exec_time = time.time()
        res = minimize(problem, algorithm, termination, seed=1, save_history=False, verbose=True)
        result = res
        exec_time = time.time() - exec_time
        self.plot_pareto(res)
        self.save_pareto_params(res)
        return result, exec_time


    @staticmethod
    def plot_pareto(res):
        plot = Scatter()
        # plot.add(problem.pareto_front(), plot_type="line", color="black", alpha=0.7)
        plot.add(res.F, facecolor="none", edgecolor="red")
        plot.show()


    @staticmethod
    def constraints(consts: list):
        consts_list = list()
        for const in consts:
            consts_list.append(eval(const))
        return consts_list


    @staticmethod
    def build_boundaries(pos_bounds):
        bounds_lower = np.array([bound[0] for bound in pos_bounds])
        bounds_upper = np.array([bound[1] for bound in pos_bounds])
        return bounds_lower, bounds_upper

    @staticmethod
    def save_pareto_params(result):
        with open(f'outputs/optimization/{date.today()}.txt', 'a') as f:
            string_to_be_saved = \
                f'{result.X}\n{result.F}\n'
            f.write(string_to_be_saved)
