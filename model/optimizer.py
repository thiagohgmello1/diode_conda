import time
import numpy as np

from model.material import Material
from model.particle import Particle
from model.topology import Topology
from simulators.monte_carlo import monte_carlo
from scipy.optimize import differential_evolution, NonlinearConstraint


class Optimizer:
    def __init__(
            self,
            pop_size: int,
            max_iter: int,
            mutation: float,
            recombination: float,
            polish: bool,
            objectives: dict,
            constraints: list,
            bounds: dict,
            geo_mask: list,
            cur_segments: list,
            material: Material,
            particle_model: Particle,
            convergence: dict,
            scale: float
    ):
        self.pop_size = pop_size
        self.max_iter = max_iter
        self.mutation = mutation
        self.material = material
        self.scale = scale
        self.geo_mask = geo_mask
        self.cur_segments = cur_segments
        self.particle_m = particle_model
        self.convergence = convergence
        self.recombination = recombination
        self.polish = polish
        self.objectives = objectives
        self.obj_func = self.choose_objective_func()
        self.consts = self.constraints(constraints)
        self.boundaries = self.build_boundaries(bounds)


    def build_geometry(self, dimensions):
        n_dim = len(dimensions)
        geo_mask = self.geo_mask
        for i in range(n_dim):
            geo_mask = [
                list(map(lambda x: dimensions[i] if x == f'{i}' else x, points_list)) for points_list in geo_mask
            ]
        return [geo_mask]


    def zero_voltage_imp(self, dimensions):
        current, voltage = self.run_specific_points(dimensions)
        impedance = 1 / np.abs(np.gradient(np.array(current), voltage))
        return impedance[1]


    def zero_bias_responsivity(self, dimensions):
        current, voltage = self.run_specific_points(dimensions)
        current_first_derivative = np.abs(np.gradient(np.array(current), voltage))
        current_second_derivative = np.abs(np.gradient(current_first_derivative, voltage))
        responsivity = current_second_derivative[1] / (2 * current_first_derivative[1])
        return 1 / np.abs(responsivity)


    def asymmetry(self, dimensions):
        pass


    def run_specific_points(self, dimensions):
        dimensions = self.integer_params(dimensions)
        topology_points = self.build_geometry(dimensions)
        topology = Topology.from_points(topology_points, self.scale, tuple(self.cur_segments))
        voltage = list()
        current = list()
        voltage_range = self.objectives['voltage_range']
        for volt in voltage_range:
            e_field, simulation_current, time_steps_count, collisions_count = \
                monte_carlo(volt, topology, self.material, self.particle_m, **self.convergence, plot_current=False)
            voltage.append(volt)
            current.append(simulation_current)
        return current, voltage


    def choose_objective_func(self):
        if self.objectives["method"] == "ZBI":
            return self.zero_voltage_imp
        elif self.objectives["method"] == "ZBR":
            return self.zero_bias_responsivity
        else:
            return self.asymmetry


    def optimize(self):
        exec_time = time.time()
        result = differential_evolution(
            self.obj_func,
            self.boundaries,
            maxiter=self.max_iter,
            popsize=self.pop_size,
            polish=self.polish,
            mutation=self.mutation,
            recombination=self.recombination,
            disp=True
        )
        exec_time = time.time() - exec_time
        return result, exec_time


    @staticmethod
    def constraints(consts: list):
        consts_list = []
        for const in consts:
            nlc = NonlinearConstraint(lambda x: eval(const[0]), eval(const[1]), eval(const[2]))
            consts_list.append(nlc)
        return consts_list


    @staticmethod
    def build_boundaries(pos_bounds):
        return [tuple(bound) for bound in pos_bounds]


    @staticmethod
    def integer_params(params):
        return [round(param) for param in params]

