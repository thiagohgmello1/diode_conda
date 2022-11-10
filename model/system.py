import numpy as np

from skgeom import Vector2, Point2, Segment2

from skgeom.draw import draw
from model.particle import Particle
from model.topology import Topology
from model.material import Material
from utils.probabilistic_operations import decision, random_number
from utils.complementary_operations import vec_to_point, calc_normal, create_segments


TEST = True


class System:
    def __init__(
            self,
            particles: list[Particle],
            topology: Topology,
            material: Material,
            electric_field: Vector2,
            delta_t: float = 1, # remover
            max_collisions: float = 1e6,
            max_time_steps: float = 1e4,
    ):
        self.particles = particles
        self.topology = topology
        self.material = material
        self.e_field = electric_field
        self.time = delta_t
        self.max_collisions = max_collisions
        self.max_time_steps = max_time_steps


    def set_init_positions(self):
        for particle in self.particles:
            particle.set_init_position(self.topology.bbox)
            init_pos = vec_to_point(particle.position)
            while not self.topology.contains(init_pos):
                particle.set_init_position(self.topology.bbox)
                init_pos = vec_to_point(particle.position)


    def relaxation_event(self, delta_t: float, collision_segment: Segment2):
        relax_time = self.material.relax_time
        relax_probability = 1 - np.exp(-delta_t / relax_time)
        relaxation = decision(relax_probability)
        if relaxation:
            delta_t = random_number(0, delta_t)
            collision_segment = None
        return relaxation, delta_t, collision_segment


    def simulate_drude(self, particle: Particle):
        collisions = 0
        time_steps = 0
        remaining_time = self.time
        condition = False
        while not condition:
            traveled_path = particle.calc_position(remaining_time)
            intersection_points = self.topology.intersection_points(traveled_path, particle.position)
            if len(intersection_points) > 0:
                collisions += 1
            lowest_time_to_collision, lowest_collision_segment = self.calc_closer_intersection(
                remaining_time, intersection_points, particle, traveled_path
            )
            particle_pos = Point2(particle.position.x(), particle.position.y())
            relaxation, lowest_time_to_collision, lowest_collision_segment = self.relaxation_event(
                lowest_time_to_collision, lowest_collision_segment
            )
            if TEST:
                particle.positions.append(particle_pos)
            segment_normal_vec = calc_normal(lowest_collision_segment, particle_pos)
            particle.move(
                segment_normal_vec, lowest_time_to_collision, self.e_field, self.material, relaxation
            )
            remaining_time -= lowest_time_to_collision
            time_steps += 1
            condition = \
                np.isclose(remaining_time, 0) or collisions > self.max_collisions or time_steps > self.max_time_steps


    @staticmethod
    def calc_closer_intersection(
            remaining_time,
            intersection_points: list[Point2],
            particle: Particle,
            traveled_path: Segment2
    ):
        lowest_time_to_collision = remaining_time
        lowest_collision_segment = None
        for intersection_point, collision_element in intersection_points:
            time_to_collision = particle.time_to_collision(intersection_point, traveled_path)
            if time_to_collision < lowest_time_to_collision:
                lowest_time_to_collision = time_to_collision
                lowest_collision_segment = collision_element
        return lowest_time_to_collision, lowest_collision_segment


if __name__ == '__main__':
    mat = Material(1, 1, 1)
    particles_list = [Particle(1, 1, 1, 10)]
    pol = Topology.from_file('../tests/test2.svg', 1e-9)
    e_field = Vector2(1, 0)
    system = System(particles_list, pol, mat, e_field, 100)
    system.set_init_positions()
    system.simulate_drude(particles_list[0])
    draw(system.topology.topologies)
    if TEST:
        segments = create_segments(particles_list[0].positions)
        draw(segments)
    print('ei')