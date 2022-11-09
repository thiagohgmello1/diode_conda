import numpy as np

from skgeom import Vector2, Point2, Segment2

from skgeom.draw import draw
from model.particle import Particle
from model.topology import Topology
from model.material import Material
from utils.complementary_operations import vec_to_point, calc_normal, create_segments
from utils.probabilistic_operations import decision, random_number


TEST = True


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


    def relaxation_event(self, time_interval: float, collision_segment: Segment2):
        relax_time = self.topology.material.relax_time
        relax_probability = 1 - np.exp(-time_interval / relax_time)
        relaxation = decision(relax_probability)
        if relaxation:
            time_interval = random_number(0, time_interval)
            collision_segment = None
        return relaxation, time_interval, collision_segment


    def simulate(self, particle: Particle):
        particle.calc_velocity(self.e_field, self.topology.material)
        while not np.isclose(self.time, 0):
            traveled_path = particle.calc_next_pos(self.time)
            intersection_points = self.topology.intersection_points(traveled_path, particle.position)
            lowest_time_to_collision, lowest_collision_segment = self.calc_closer_intersection(
                intersection_points, particle, traveled_path
            )
            particle_pos = Point2(particle.position.x(), particle.position.y())
            relaxation, lowest_time_to_collision, lowest_collision_segment = self.relaxation_event(
                lowest_time_to_collision, lowest_collision_segment
            )
            if TEST:
                particle.positions.append(particle_pos)
            segment_normal_vec = calc_normal(lowest_collision_segment, particle_pos)
            particle.move(
                segment_normal_vec, lowest_time_to_collision, self.e_field, self.topology.material, relaxation
            )
            self.time -= lowest_time_to_collision


    def calc_closer_intersection(self, intersection_points: list[Point2], particle: Particle, traveled_path: Segment2):
        lowest_time_to_collision = self.time
        lowest_collision_segment = None
        for intersection_point, collision_element in intersection_points:
            time_to_collision = particle.time_to_collision(intersection_point, traveled_path)
            if time_to_collision < lowest_time_to_collision:
                lowest_time_to_collision = time_to_collision
                lowest_collision_segment = collision_element
        return lowest_time_to_collision, lowest_collision_segment


if __name__ == '__main__':
    particles_list = [Particle(1, 1, 1, 10)]
    mat = Material(1, 1, 1)
    pol = Topology.from_file('../tests/test2.svg', mat)
    e_field = Vector2(1, 0)
    system = System(particles_list, pol, e_field, 100)
    system.set_init_positions()
    system.simulate(particles_list[0])
    draw(system.topology.topologies)
    if TEST:
        segments = create_segments(particles_list[0].positions)
        draw(segments)
    print('ei')
