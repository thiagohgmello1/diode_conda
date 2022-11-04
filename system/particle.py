from utils.random_gen import random_vec
from utils.vec_operations import mirror
from skgeom import Bbox2, Vector2
from material import Material

import numpy as np


class Particle:
    def __init__(self, charge: float, mass: float, size: float, fermi_velocity: float, position=None):
        self.charge = charge
        self.mass = mass
        self.size = size
        self.fermi_velocity = random_vec() * fermi_velocity
        self.velocity = None
        self.position = position
        self.path = None


    def calc_velocity(self, electric_field: np.array, material: Material):
        avg_vel = -1 * (electric_field * self.charge * material.relax_time) / (self.mass * material.mass_scale)
        self.velocity = self.fermi_velocity + avg_vel


    def set_init_position(self, bbox: Bbox2):
        min_range = (bbox.xmin(), bbox.ymin())
        max_range = (bbox.xmax(), bbox.ymax())
        self.position = random_vec(min_value=min_range, max_value=max_range, is_normalized=False)


    def calc_next_pos(self, time: float):
        return self.position + self.velocity * time


    def move(self, normal_vec: np.array, time: float, electric_field: np.array, material: Material):
        self.position = self.position + self.velocity * time
        self.velocity = mirror(self.velocity, normal_vec)
        self.calc_velocity(electric_field, material)


if __name__ == '__main__':
    particle = Particle(1, 1, 1, 10)
    mat = Material(1, 1)
    e_field = np.array([1, 0])
    particle.calc_velocity(e_field, mat)
    box = Bbox2(1, 2, 3, 4)
    particle.set_init_position(box)
    print('ei')
