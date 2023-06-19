import numpy as np
from skgeom import Vector2


class ElectricField:
    def __init__(self):
        self.max_mag = 0
        self.field_dir = np.array([1, 0])
        self.vector = lambda vec: self.max_mag


    def calc_electric_field(self, voltage: float, method: str, device_len: float):
        if method == 'consistent':
            # Implement self-consistent field calculation based on Poisson equation
            pass
        else:
            e_field_aux = voltage * self.field_dir
            self.max_mag = Vector2(*e_field_aux) / device_len

