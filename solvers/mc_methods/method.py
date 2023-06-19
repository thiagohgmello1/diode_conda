from math import floor, log10
from abc import ABC, abstractmethod
from models.particle import Particle

SIGNIFICANT_DIGITS = 4


class Method(ABC):
    def __init__(self, system):
        self.system = system
        self.voltage = None
        self.significant_digits_time = -(int(floor(log10(self.system.material.relax_time)))) + SIGNIFICANT_DIGITS

    @abstractmethod
    def simulate(self, particle: Particle):
        pass


    @abstractmethod
    def calc_current(self):
        pass
