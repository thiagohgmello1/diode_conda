import numpy as np

from datetime import datetime
from model.material import Material
from model.particle import Particle
from model.topology import Topology
from simulators.monte_carlo import monte_carlo


class Optimizer:
    def __init__(
            self,
            pop_size: int,
            max_iter: int,
            objectives: dict,
            geo_mask: list,
            cur_segments: list,
            material: Material,
            particle_model: Particle,
            convergence: dict,
            scale: float
    ):
        self.pop_size = pop_size
        self.max_iter = max_iter
        self.material = material
        self.scale = scale
        self.result = tuple()
        self.geo_mask = geo_mask
        self.cur_segments = cur_segments
        self.particle_m = particle_model
        self.convergence = convergence
        self.objectives = objectives
        self.derivative_tech = self.def_derivative_technique()


    def def_derivative_technique(self):
        if self.objectives['derivative'] == 'fit':
            return self.poly_fit_derivatives_zero_bias
        else:
            return self.numerical_derivatives_zero_bias


    def build_geometry(self, dimensions):
        n_dim = len(dimensions)
        geo_mask = self.geo_mask
        for i in range(n_dim):
            geo_mask = [
                list(
                    map(
                        lambda x: x.replace(f'x{i}', str(dimensions[i]))
                        if (isinstance(x, str) and f'x{i}' in x) else x, points_list
                    )
                ) for points_list in geo_mask
            ]
        geo_mask = list([*map(lambda x: eval(x) if isinstance(x, str) else x, points_list)] for points_list in geo_mask)
        return [geo_mask]


    def zero_voltage_imp(self, dimensions):
        if ("method" in self.objectives.keys()) or (self.objectives["methods"][0] == "ZBI"):
            self.result = self.run_specific_points(dimensions)
        current_first_derivative, _ = self.derivative_tech(self.result[0], self.result[1])
        impedance = 1 / current_first_derivative
        return impedance


    def zero_bias_responsivity(self, dimensions):
        if ("method" in self.objectives.keys()) or (self.objectives["methods"][0] == "ZBR"):
            self.result = self.run_specific_points(dimensions)
        current_first_derivative, current_second_derivative = self.derivative_tech(self.result[0], self.result[1])
        responsivity = current_second_derivative / (2 * current_first_derivative)
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


    def poly_fit_derivatives_zero_bias(self, current, voltage):
        iv_f = np.poly1d(np.polyfit(voltage, current, self.objectives['poly_order']))
        current_first_derivative = np.polyder(iv_f, 1)
        current_second_derivative = np.polyder(iv_f, 2)
        return current_first_derivative(0), current_second_derivative(0)


    @staticmethod
    def integer_params(params):
        return [round(param) for param in params]


    @staticmethod
    def save_current_iter(x, convergence):
        now = datetime.now()
        date_string = now.strftime("%d/%m/%Y %H:%M:%S")
        with open(f'outputs/optimization/opt_intermediates.csv', 'a') as f:
            string_to_be_saved = \
                f'{date_string};{x};{convergence}\n'
            f.write(string_to_be_saved)


    @staticmethod
    def numerical_derivatives_zero_bias(current, voltage):
        zero_index = voltage.index(0)
        current_first_derivative = np.abs(np.gradient(np.array(current), voltage))
        current_second_derivative = np.abs(np.gradient(current_first_derivative, voltage))
        return current_first_derivative[zero_index], current_second_derivative[zero_index]
