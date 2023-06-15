import matplotlib
import numpy as np
import multiprocessing

from math import log10, floor
from model.particle import Particle
from model.topology import Topology
from model.material import Material
from skgeom import Vector2, Point2, Segment2
from scipy.constants import elementary_charge
from utils.post_processing import progress_bar, plot_stable_current
from utils.complementary_operations import vec_to_point, point_to_vec, calc_normal


TEST = False
followed_particle_id = 1
TIME_PRECISION = 0.99
SIGNIFICANT_DIGITS = 4
matplotlib.use('TkAgg')


class System:
    def __init__(
            self,
            particle: Particle,
            topology: Topology,
            material: Material,
            electric_field: Vector2,
            number_of_particles: int = None,
            max_collisions: float = np.inf,
            max_time_steps: float = np.inf
    ):
        """
        Create system to be simulated (topology + particles + materials + etc.)

        :param particle: particle model to be simulated
        :param topology: desired topology
        :param material: material
        :param electric_field: defined or calculated electric field created by applied topology voltage [V/m]
        :param number_of_particles: number of particles (defined by the number of cores if not defined)
        :param max_collisions: defined maximum accepted collisions. Stop criteria
        :param max_time_steps: defined maximum time steps. Stop criteria [s]
        """
        self.currents = list()

        self.particles = self.create_particles(particle, number_of_particles)
        self.topology = topology
        self.material = material
        self.e_field = electric_field
        self.relax_time = self.material.relax_time

        self.total_macro_particles = len(self.particles) * particle.density
        self.particles_counter = 0

        self.simulated_time = 0

        self.max_collisions = max_collisions
        self.collisions_count = 0

        self.max_time_steps = max_time_steps
        self.time_steps_count = 0

        self.significant_digits_time = -(int(floor(log10(self.material.relax_time)))) + SIGNIFICANT_DIGITS


    @staticmethod
    def create_particles(particle_model, number_of_particles) -> list:
        """
        Create particles according specified number or according to CPU cores

        :param particle_model: model of desired particle
        :param number_of_particles: number of particles to be created (if specified)
        :return: list of particles
        """
        particles = list()
        if not number_of_particles:
            number_of_particles = multiprocessing.cpu_count()
        for _ in range(number_of_particles):
            particles.append(particle_model.__deepcopy__())
        return particles


    def set_particle_parameters(self):
        """
        Set particle initial parameters

        :return: None
        """
        for particle in self.particles:
            particle.set_init_position(self.topology.bbox)
            init_pos = vec_to_point(particle.position)
            _, lowest_dist = self.topology.get_closer_segment(init_pos)
            while (not self.topology.contains(init_pos)) or (lowest_dist < self.topology.scale / 20):
                particle.set_init_position(self.topology.bbox)
                init_pos = vec_to_point(particle.position)
                _, lowest_dist = self.topology.get_closer_segment(init_pos)
            if TEST and particle.id == followed_particle_id:
                particle.positions.append([init_pos])


    def simulate(self, model, voltage: list):
        """
        Simulate complete system. Simulation can be executed in GPU or single-multi core CPU

        :param voltage: simulated voltage
        :param model: desired model to simulate. Must be a method callback (ex.: simulate_drude method)
        :return: None
        """
        self.set_particle_parameters()
        while not self._stop_conditions():
            self.simulated_time += self.relax_time
            self.time_steps_count += 1
            for particle in self.particles:
                particle.set_velocity()
                model(particle)
                if self._stop_conditions():
                    break
            self.currents.append(self.cal_current())
            progress_bar(self.collisions_count, self.max_collisions)
        # plot_stable_current(self.currents, voltage)
        print('\n')


    def _stop_conditions(self) -> bool:
        """
        Calculate if any stop condition was met

        :return: stop condition
        """
        time_steps_condition = self.time_steps_count > self.max_time_steps
        collisions_condition = self.collisions_count > self.max_collisions

        return time_steps_condition or collisions_condition


    def simulate_drude(self, particle: Particle):
        """
        Simulate Drude event

        :param particle: macroparticle to be simulated
        :return: None
        """
        drift_velocity = particle.calc_drift_velocity(self.relax_time, self.e_field, self.material.mobility, 'relax')
        particle.velocity += drift_velocity

        remaining_time = self.relax_time
        stop_conditions = np.isclose(
            remaining_time, 0, atol=10 ** (-self.significant_digits_time)) or remaining_time < 0

        while not stop_conditions:
            traveled_path = particle.calc_next_position(remaining_time)
            intersection_points = self.topology.intersection_points(traveled_path)
            lowest_time_to_collision, closest_collision_segment, next_pos = self.calc_closer_intersection(
                remaining_time, particle.velocity, intersection_points, traveled_path, particle
            )
            particle_p0 = vec_to_point(particle.position)
            particle.position = next_pos

            if not self.topology.contains(vec_to_point(particle.position)):
                raise Exception('Particle is outside geometry')

            collision_normal_vec = calc_normal(closest_collision_segment, particle_p0)
            self.check_current_segment_collision(particle, closest_collision_segment, collision_normal_vec)
            remaining_time -= lowest_time_to_collision
            self.save_particle_data(particle)
            stop_conditions = np.isclose(
                remaining_time, 0, atol=10 ** (-self.significant_digits_time)) or remaining_time < 0


    def check_current_segment_collision(self, particle, closest_collision_segment, segment_normal_vec):
        """
        Check if macroparticle collide with current computation segment

        :param particle: travelled particle
        :param closest_collision_segment: collided segment
        :param segment_normal_vec: collided segment normal vector
        :return: None
        """
        if segment_normal_vec:
            self.collisions_count += 1
            current_collision, element = self.particle_computation(closest_collision_segment)
            if current_collision:
                self.teleport_particle(particle, element)
            else:
                particle.mirror_particle(segment_normal_vec)


    def particle_computation(self, collided_element: Segment2) -> tuple[bool, str]:
        """
        Compute collisions in current elements

        :param collided_element: collided geometry segment
        :return: boolean indicating if there was rectification
        :return: string indicating segment group
        """
        if collided_element and collided_element in self.topology.current_computing_elements['direct']:
            self.particles_counter += self.particles[0].density
            return True, 'reverse'
        elif collided_element and collided_element in self.topology.current_computing_elements['reverse']:
            self.particles_counter -= self.particles[0].density
            return True, 'direct'
        return False, ''


    def teleport_particle(self, particle: Particle, element: str):
        """
        Teleport particle to opposite current segment
        :param particle: travelled macroparticle
        :param element: segment current group (i.e. 'direct' or 'reverse')
        :return: None
        """
        pos = self.topology.random_segment_pos(element)
        # pos = Point2(pos.x(), particle.position.y())
        particle.position = point_to_vec(pos)


    def calc_closer_intersection(
            self,
            remaining_time,
            particle_velocity: Vector2,
            intersection_points: list[Point2],
            traveled_path: Segment2,
            particle: Particle
    ) -> (float, Segment2, Vector2):
        """
        Define closer intersection between particle path and geometries boundaries

        :param remaining_time: time until scattering process
        :param particle_velocity: possible particle velocity
        :param intersection_points: list of points where collision can occur
        :param traveled_path: corresponding particle path
        :param particle: simulated particle
        :return lowest_time_to_collision: the lowest time to collision
        :return lowest_collision_segment: collided segment
        :return next_pos: particle next position
        """
        lowest_time_to_collision = remaining_time
        next_pos = None
        lowest_collision_segment = None
        for intersection_point, collision_element in intersection_points:
            time_to_collision = self._time_to_collision(particle_velocity, intersection_point, traveled_path)
            if time_to_collision < lowest_time_to_collision:
                lowest_time_to_collision = time_to_collision * TIME_PRECISION
                lowest_collision_segment = collision_element
                next_pos = intersection_point

        if not next_pos:
            next_pos = particle.position + particle.velocity * lowest_time_to_collision

        return lowest_time_to_collision, lowest_collision_segment, next_pos


    def cal_current(self):
        """
        Calculate total current

        :return: calculated current
        """
        carrier_concentration = self.material.carrier_concentration
        current = (
                          carrier_concentration * self.topology.area * elementary_charge * self.particles_counter
                  ) / (self.simulated_time * self.total_macro_particles)
        return current


    @staticmethod
    def _time_to_collision(particle_velocity: Vector2, position: Point2, path: Segment2) -> float:
        """
        Calculate time until collision happen. Consider uniform particle movement in delta_t

        :param particle_velocity: particle velocity
        :param position: particle possible next position
        :param path: estimated particle path
        :return: time until collision point
        """
        pos = Segment2(path[0], position)
        return np.sqrt(float(pos.squared_length() / particle_velocity.squared_length()))


    @staticmethod
    def save_particle_data(particle):
        """
        Save specific particle positions (if TEST is equal True)

        :param particle: desired particle to save data
        :return: None
        """
        if TEST and particle.id == followed_particle_id:
            particle_pos = Point2(particle.position.x(), particle.position.y())
            particle.positions.append(particle_pos)
