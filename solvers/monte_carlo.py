import time
import numpy as np

from solvers.solver import Solver
from scipy.constants import electron_mass
from matplotlib.ticker import EngFormatter
from solvers.mc_methods.method import Method
from solvers.drude_analytical import drude_analytical_model


class MonteCarlo(Solver):
    def __init__(self, method: Method, max_coll: int = np.inf, max_time_steps: int = np.inf, geo: str = None):
        super().__init__(method, geo)

        self.max_collisions = max_coll
        self.max_time_steps = max_time_steps
        self.eng_formatter = EngFormatter(places=4, unit='A')


    def stop_conditions(self) -> bool:
        """
        Calculate if any stop condition was met

        :return: stop condition
        """
        collisions_condition = self.system.collisions_count > self.max_collisions
        time_steps_condition = self.system.time_steps_count > self.max_time_steps

        return time_steps_condition or collisions_condition


    def simulate(self):
        for voltage in self.voltage_range:
            self.system.set_particles_parameters()
            device_len = self.system.topology.bbox.xmax() - self.system.topology.bbox.xmin()
            self.system.electric_field.calc_electric_field(voltage, 'static', device_len)
            self.method.voltage = voltage
            partial_exec_times = time.time()
            self.simulate_stop_cond()

            self.voltages.append(voltage)
            calc_current = self.method.calc_current()
            self.currents.append(calc_current)
            partial_exec_times = time.time() - partial_exec_times
            self.partial_exec_times.append(partial_exec_times)
            self.compare_sim()
            self.print_data()
            self.system.reset_performance_params()
        return self.voltages, self.currents, self.compare_currents, self.partial_exec_times


    def simulate_stop_cond(self):
        """
        iterate voltage until some stop condition is not reached
        :return:
        """
        while not self.stop_conditions():
            self.system.simulated_time += self.system.material.relax_time
            self.system.time_steps_count += 1
            self.simulate_particles()
            self.progress_bar(self.system.collisions_count, self.max_collisions)
        print('\r')


    def simulate_particles(self):
        """
        Simulate all particles

        :return:
        """
        self.system.scatter_particles()
        for particle in self.system.particles.particles_list:
            self.method.simulate(particle)
            if self.stop_conditions():
                break


    def compare_sim(self):
        if 'rectangle' in self.geo:
            drude_current = drude_analytical_model(
                e_field=self.system.electric_field,
                relax_time=self.system.material.relax_time,
                width=float(self.system.topology.bbox.ymax() - self.system.topology.bbox.ymin()),
                carrier_concentration=self.system.material.carrier_concentration,
                effective_mass=self.system.material.effective_mass * electron_mass
            )
            self.compare_currents.append(drude_current)


    def print_data(self):
        print('\n')
        print(f"Voltage: {'%s' % float('%.1g' % self.voltages[-1])}")
        print(f"Current:{self.eng_formatter.format_eng(num=self.currents[-1])}A")
        if 'rectangle' in self.geo:
            print(f'Drude current: {self.eng_formatter.format_eng(num=self.compare_currents[-1])}A')
        print(f'Time steps: {self.system.time_steps_count}')
        print(f'Collisions: {self.system.collisions_count}')
        print(f'-' * 100)
        print('\r')



