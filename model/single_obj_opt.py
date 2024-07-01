import time
import numpy as np

from model.material import Material
from model.optimizer import Optimizer
from scipy.optimize import differential_evolution, NonlinearConstraint

from model.particle import Particle


class SingleObjOpt(Optimizer):
    def __init__(self, params: dict, material: Material, particle_model: Particle, convergence: dict, scale: float):
        self.mutation = params.pop("mutation")
        self.polish = params.pop("polish")
        self.recombination = params.pop("recombination")
        self.consts = self.constraints(params.pop("constraints"))
        self.boundaries = self.build_boundaries(params.pop("bounds"))
        super().__init__(
            **params, material=material, particle_model=particle_model, convergence=convergence, scale=scale
        )
        self.obj_funcs = self.choose_objective_func()

    def optimize(self):
        exec_time = time.time()
        result = differential_evolution(
            self.obj_funcs,
            self.boundaries,
            maxiter=self.max_iter,
            popsize=self.pop_size,
            polish=self.polish,
            mutation=self.mutation,
            recombination=self.recombination,
            disp=True,
            constraints=self.consts,
            callback=self.save_current_iter
        )
        exec_time = time.time() - exec_time
        return result, exec_time


    def choose_objective_func(self):
        if self.objectives["method"] == "ZBI":
            return self.zero_voltage_imp
        elif self.objectives["method"] == "ZBR":
            return self.zero_bias_responsivity
        else:
            return self.asymmetry


    @staticmethod
    def constraints(consts: list):
        consts_set = set()
        for const in consts:
            nlc = NonlinearConstraint(eval(const[0]), eval(const[1]), eval(const[2]))
            consts_set.add(nlc)
        return consts_set


    @staticmethod
    def build_boundaries(pos_bounds):
        return [tuple(bound) for bound in pos_bounds]
