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
            delta_t: float,
            max_time_simulation: float = np.inf,
            max_collisions: float = np.inf,
            max_time_steps: float = np.inf,
    ):
        """
        Create system to be simulated (topology + particles + materials + etc.)

        :param particles: particles to be simulated
        :param topology: desired topology
        :param material: material to topology
        :param electric_field: defined or calculated electric field created by applied voltage in topology
        :param delta_t: time step
        :param max_collisions: defined maximum accepted collisions. Stop criteria
        :param max_time_steps: defined maximum time steps. Stop criteria
        """
        self.particles = particles
        self.topology = topology
        self.material = material
        self.e_field = electric_field
        self.delta_t = delta_t

        self.collisions = 0
        self.time_steps = 0

        self.max_time_simulation = max_time_simulation
        self.max_collisions = max_collisions
        self.max_time_steps = max_time_steps


    def set_particles_parameters(self):
        """
        Set particles initial parameters

        :return: None
        """
        for particle in self.particles:
            particle.set_fermi_velocity()
            particle.set_init_position(self.topology.bbox)
            init_pos = vec_to_point(particle.position)

            while not self.topology.contains(init_pos):
                particle.set_init_position(self.topology.bbox)
                init_pos = vec_to_point(particle.position)


    def relaxation_event(self, delta_t: float, collision_segment: Segment2) -> (bool, float, Segment2):
        """
        Define if relaxation event occur in time interval

        :param delta_t: time interval
        :param collision_segment: segment where particle will collide
        :return relaxation: boolean indicating if there is or there is not relaxation event
        :return delta_t: time interval
        :return collision_segment: segment where particle will collide if there is no relaxation event
        """
        relax_probability = 1 - np.exp(-delta_t / self.material.relax_time)
        relaxation = decision(relax_probability)

        if relaxation:
            delta_t = random_number(0, delta_t)
            collision_segment = None

        return relaxation, delta_t, collision_segment


    def simulate_drude(self, particle: Particle):
        """
        Simulate Drude event. Simulation can be executed in GPU or single-multithreading CPU

        :param particle: particle to be simulated
        :return: None
        """
        simulated_time = 0
        condition = self.stop_condition(simulated_time)

        while not condition:
            particle.calc_velocity(self.e_field, self.delta_t)
            traveled_path = particle.calc_next_position(self.delta_t)
            topology_contains_next_point = self.topology.contains(traveled_path[1])

            if not topology_contains_next_point:
                intersection_points = self.topology.intersection_points(traveled_path)
                self.collisions += 1
            else:
                intersection_points = list()

            lowest_time_to_collision, lowest_collision_segment = self.calc_closer_intersection(
                intersection_points, particle, traveled_path
            )
            particle_pos = Point2(particle.position.x(), particle.position.y())
            relaxation, lowest_time_to_collision, lowest_collision_segment = self.relaxation_event(
                lowest_time_to_collision, lowest_collision_segment
            )
            segment_normal_vec = calc_normal(lowest_collision_segment, particle_pos)
            particle.move(segment_normal_vec, lowest_time_to_collision, self.e_field, relaxation)

            if relaxation:
                particle.calc_velocity(self.e_field, lowest_time_to_collision)
                if not topology_contains_next_point:
                    self.collisions -= 1

            simulated_time += lowest_time_to_collision
            self.time_steps += 1
            condition = self.stop_condition(simulated_time)
            if TEST:
                particle.positions.append(particle_pos)


    def stop_condition(self, simulated_time) -> bool:
        """
        Calculate if any stop condition was met

        :param simulated_time: total simulated time for specific particle
        :return: stop condition
        """
        time_steps_condition = self.time_steps > self.max_time_steps
        collisions_condition = self.collisions > self.max_collisions
        time_condition = np.isclose(simulated_time, self.max_time_simulation)

        return time_steps_condition or collisions_condition or time_condition


    def calc_closer_intersection(
            self,
            intersection_points: list[Point2],
            particle: Particle,
            traveled_path: Segment2
    ) -> (float, Segment2):
        """
        Define closer intersection between particle path and geometries boundaries

        :param intersection_points: list of points where collision can occur
        :param particle: respective moving particle
        :param traveled_path: corresponding particle path
        :return lowest_time_to_collision: lowest time to collision
        :return lowest_collision_segment: collided segment
        """
        lowest_time_to_collision = self.delta_t
        lowest_collision_segment = None
        for intersection_point, collision_element in intersection_points:
            time_to_collision = particle.time_to_collision(intersection_point, traveled_path)
            if time_to_collision < lowest_time_to_collision:
                lowest_time_to_collision = time_to_collision
                lowest_collision_segment = collision_element

        return lowest_time_to_collision, lowest_collision_segment


if __name__ == '__main__':
    mat = Material(1, 1, 1)
    particles_list = [Particle(1, 1, 1)]
    if TEST:
        particles_list[0].charge = 1
        particles_list[0].mass = 1
        mat.carrier_concentration = 1
        mat.effective_mass = 1
    pol = Topology.from_file('../tests/test3.svg', 1)
    e_field = Vector2(5, 0)
    system = System(particles_list, pol, mat, e_field, 10)
    system.set_particles_parameters()
    system.simulate_drude(particles_list[0])
    draw(system.topology.topologies)
    if TEST:
        segments = create_segments(particles_list[0].positions)
        draw(segments)
    print('ei')
