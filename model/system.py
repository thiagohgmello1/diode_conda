import numpy as np

from skgeom import Vector2, Point2
from utils.complementary_operations import vec_to_point, calc_normal

from skgeom.draw import draw
from model.particle import Particle
from model.topology import Topology
from model.material import Material


class System:
    def __init__(self, particles: list[Particle], topology: Topology, electric_field: Vector2, time: float):
        self.particles = particles
        self.topology = topology
        self.e_field = electric_field
        self.time = time


    def set_init_positions(self):
        for particle in self.particles:
            particle.set_init_position(self.topology.bbox)
            init_pos = vec_to_point(particle.position)
            while not self.topology.contains(init_pos):
                particle.set_init_position(self.topology.bbox)
                init_pos = vec_to_point(particle.position)


    def simulate(self, particle):
        remaining_time = self.time
        particle.calc_velocity(self.e_field, self.topology.material)
        while not np.isclose(remaining_time, 0):
            traveled_path = particle.calc_next_pos(remaining_time)
            intersection_points = self.topology.intersection_points(traveled_path, particle.position)
            lowest_time_to_collision = remaining_time
            lowest_collision_segment = None
            for intersection_point, collision_element in intersection_points:
                time_to_collision = particle.time_to_collision(intersection_point, traveled_path)
                if time_to_collision < lowest_time_to_collision:
                    lowest_time_to_collision = time_to_collision
                    lowest_collision_segment = collision_element
            remaining_time -= lowest_time_to_collision
            particle_pos = Point2(particle.position.x(), particle.position.y())
            particle.positions.append(particle_pos)
            segment_normal_vec = calc_normal(lowest_collision_segment, particle_pos)
            particle.move(segment_normal_vec, lowest_time_to_collision, self.e_field, self.topology.material)


if __name__ == '__main__':
    particles = [Particle(1, 1, 1, 10)]
    mat = Material(1, 1)
    pol = Topology.from_file('../tests/test2.svg', mat)
    e_field = Vector2(1, 0)
    system = System(particles, pol, e_field, 100)
    system.set_init_positions()
    system.simulate(particles[0])
    print('ei')
