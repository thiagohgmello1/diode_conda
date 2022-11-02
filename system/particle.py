from utils.random_gen import random_vec


class Particle:
    def __init__(self, charge, mass, size, velocity=None, position=None):
        self.charge = charge
        self.mass = mass
        self.size = size
        self.velocity = velocity
        self.position = position


    def set_velocity(self):
        self.velocity = random_vec()


