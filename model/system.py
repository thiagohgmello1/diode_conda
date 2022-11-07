from skgeom import Vector2
from utils.vec_operations import vec_to_point

from skgeom.draw import draw
from model.particle import Particle
from model.topology import Topology
from model.material import Material


class System:
    def __init__(self, particles: list[Particle], topology: Topology, time: float):
        self.particles = particles
        self.topology = topology
        self.time = time


    def set_init_positions(self):
        for particle in self.particles:
            particle.set_init_position(self.topology.bbox)
            init_pos = vec_to_point(particle.position)
            while not self.topology.contains(init_pos):
                particle.set_init_position(self.topology.bbox)
                init_pos = vec_to_point(particle.position)


    def simulate(self):
        pass


if __name__ == '__main__':
    pol = Topology.from_file('../tests/test2.svg')
    particles = [Particle(1, 1, 1, 10)]
    mat = Material(1, 1)
    e_field = Vector2(1, 0)
    system = System(particles, pol, 1)
    system.set_init_positions()
    print('ei')
