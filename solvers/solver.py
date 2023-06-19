import numpy as np

from abc import ABC, abstractmethod
from solvers.mc_methods.method import Method


class Solver(ABC):
    def __init__(self, method: Method, geo):
        self.partial_exec_times = list()
        self.currents = list()
        self.compare_currents = list()
        self.voltages = list()

        self.geo = geo
        self.method = method
        self.voltage_range = None
        self.system = self.method.system


    @abstractmethod
    def simulate(self):
        pass


    @abstractmethod
    def stop_conditions(self):
        pass


    def create_voltage_range(self, v_min: float, v_max: float, num_points: int):
        if num_points == 1:
            self.voltage_range = [v_max]
        self.voltage_range = np.linspace(v_min, v_max, num=num_points)


    @staticmethod
    def progress_bar(progress, total):
        """
        Print progress bar

        :param progress: evolved situation
        :param total: expect max situation
        :return: None
        """
        percent = 100 * (progress / total)
        bar = '█' * int(percent) + '-' * (100 - int(percent))
        print(f'\r|{bar}| {percent:.2f}%', end='\r')
