import matplotlib
import numpy as np

from math import log10, floor
from model.particle import Particle
from model.topology import Topology
from model.material import Material
from skgeom import Vector2, Point2, Segment2
from scipy.constants import elementary_charge
from utils.post_processing import progress_bar, plot_stable_current
from utils.complementary_operations import vec_to_point, point_to_vec, calc_normal, norm, calc_versor


TEST = False
followed_particle_id = 1
TIME_PRECISION = 0.99
DIST_PRECISION = 1
SIGNIFICANT_DIGITS = 4
MAX_LOOP = 200
MAX_RECONFIG = 200
matplotlib.use('TkAgg')


class System:
    def __init__(
            self,
            particle: Particle,
            topology: Topology,
            material: Material,
            electric_field: Vector2,
            check_condition: str,
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
        :param check_condition: define the method to check the scattering behaviour ('time' or 'distance')
        :param number_of_particles: number of particles (defined by the number of cores if not defined)
        :param max_collisions: defined maximum accepted collisions. Stop criteria
        :param max_time_steps: defined maximum time steps. Stop criteria [s]
        """
        self.currents = list()
        self.check_condition = check_condition

        self.particle = particle
        self.topology = topology
        self.material = material
        self.e_field = electric_field
        self.relax_time = self.material.relax_time

        self.total_macro_particles = particle.density
        self.particles_counter = 0
        self.collisions_count = 0
        self.time_steps_count = 0
        self.simulated_time = 0

        self.max_collisions = max_collisions
        self.collisions_count = 0
        self.max_time_steps = max_time_steps
        self.time_steps_count = 0

        self.significant_digits_time = -(int(floor(log10(self.material.relax_time)))) + SIGNIFICANT_DIGITS
        self.significant_digits_dist = -(int(floor(log10(self.material.mean_free_path)))) + SIGNIFICANT_DIGITS


    def set_particle_parameters(self):
        """
        Set particle initial parameters

        :return: None
        """
        self.particle.set_init_position(self.topology.bbox)
        init_pos = vec_to_point(self.particle.position)
        _, lowest_dist = self.topology.get_closer_segment(init_pos)
        while (not self.topology.contains(init_pos)) or (lowest_dist < self.topology.scale / 20):
            self.particle.set_init_position(self.topology.bbox)
            init_pos = vec_to_point(self.particle.position)
            _, lowest_dist = self.topology.get_closer_segment(init_pos)
        if TEST and self.particle.id == followed_particle_id:
            self.particle.positions.append(init_pos)


    def simulate(self, model, voltage: list, plot_current: bool = True):
        """
        Simulate complete system

        :param voltage: simulated voltage
        :param model: desired model to simulate. Must be a method callback (ex.: simulate_drude method)
        :param plot_current: define if stable current will be plotted
        :return: None
        """
        reconfig_count = 0
        self.set_particle_parameters()
        while (not self._stop_conditions()) and (reconfig_count <= MAX_RECONFIG):
            self.time_steps_count += 1
            self.particle.set_velocity()
            traveled_time, loop_condition = model()
            if loop_condition:
                self.set_particle_parameters()
                reconfig_count += 1
            else:
                self.simulated_time += traveled_time
                if self._stop_conditions():
                    break
                self.currents.append(self.cal_current())
                progress_bar(self.collisions_count, self.max_collisions)
        if reconfig_count > MAX_RECONFIG:
            raise Exception('Max reconfiguration')
        if plot_current:
            plot_stable_current(self.currents, voltage)
        print('\n')


    def _stop_conditions(self) -> bool:
        """
        Calculate if any stop condition was met

        :return: stop condition
        """
        time_steps_condition = self.time_steps_count > self.max_time_steps
        collisions_condition = self.collisions_count > self.max_collisions

        return time_steps_condition or collisions_condition


    def simulate_drude(self):
        """
        Simulate Drude event

        :return: None
        """
        self.particle.calc_drift_velocity(self.relax_time, self.e_field, self.material.mobility)
        count_loop = 0

        traveled_time = 0
        remaining_time = self.relax_time
        remaining_dist = self.material.mean_free_path
        stop_conditions, loop_cond = self._calc_stop_conditions(remaining_time, remaining_dist, count_loop)

        while not (stop_conditions or loop_cond):
            traveled_path = self.calc_traveled_path(remaining_time, remaining_dist)
            intersection_points = self.topology.intersection_points(traveled_path)
            lowest_time_to_collision, lowest_dist_to_collision, closest_collision_segment, next_pos = \
                self.calc_closer_intersection(remaining_time, remaining_dist, intersection_points, traveled_path)
            particle_p0 = vec_to_point(self.particle.position)
            self.particle.position = next_pos

            if not self.topology.contains(vec_to_point(self.particle.position)):
                raise Exception('Particle is outside geometry')

            collision_normal_vec = calc_normal(closest_collision_segment, particle_p0)
            self.check_current_segment_collision(closest_collision_segment, collision_normal_vec)
            traveled_time += lowest_time_to_collision
            remaining_time -= lowest_time_to_collision
            remaining_dist -= lowest_dist_to_collision
            count_loop += 1
            stop_conditions, loop_cond = self._calc_stop_conditions(remaining_time, remaining_dist, count_loop)
            self.save_particle_data()

        return traveled_time, loop_cond


    def calc_traveled_path(self, delta_t: float, delta_s: float) -> Segment2:
        next_pos = self.particle.calc_next_position(delta_t, delta_s, self.check_condition)
        p_0 = vec_to_point(self.particle.position)
        p_1 = vec_to_point(next_pos)
        travel_path = Segment2(p_0, p_1)

        return travel_path


    def check_current_segment_collision(self, closest_collision_segment, segment_normal_vec):
        """
        Check if macroparticle collide with current computation segment

        :param closest_collision_segment: collided segment
        :param segment_normal_vec: collided segment normal vector
        :return: None
        """
        if segment_normal_vec:
            self.collisions_count += 1
            current_collision, element = self.particle_computation(closest_collision_segment)
            if current_collision:
                self.teleport_particle(element)
            else:
                self.particle.mirror_particle(segment_normal_vec)


    def particle_computation(self, collided_element: Segment2) -> tuple[bool, str]:
        """
        Compute collisions in current elements

        :param collided_element: collided geometry segment
        :return: boolean indicating if there was rectification
        :return: string indicating segment group
        """
        if collided_element and collided_element in self.topology.current_computing_elements['direct']:
            self.particles_counter += self.particle.density
            return True, 'reverse'
        elif collided_element and collided_element in self.topology.current_computing_elements['reverse']:
            self.particles_counter -= self.particle.density
            return True, 'direct'
        return False, ''


    def teleport_particle(self, element: str):
        """
        Teleport particle to opposite current segment

        :param element: segment current group (i.e. 'direct' or 'reverse')
        :return: None
        """
        pos = self.topology.random_segment_pos(element)
        # pos = self.topology.specific_segment_pos(element, self.particle.position)
        # pos = Point2(pos.x(), particle.position.y())
        self.particle.position = point_to_vec(pos)


    def calc_closer_intersection(
            self,
            remaining_time,
            remaining_dist,
            intersection_points: list[Point2],
            traveled_path: Segment2
    ) -> (float, float, Segment2, Vector2):
        """
        Calculate closer intersection between particle path and geometries boundaries according to stop criterion

        :param remaining_time: time until scattering process
        :param remaining_dist: distance until scattering process
        :param intersection_points: list of points where collision can occur
        :param traveled_path: corresponding particle path
        :return lowest_time_to_collision: the lowest time to collision
        :return lowest_collision_segment: collided segment
        :return next_pos: particle next position
        """

        if self.check_condition == "time":
            lowest_time_to_collision, lowest_collision_segment, next_pos = \
                self.time_intersection(remaining_time, intersection_points, traveled_path)
            lowest_dist_to_collision = norm(calc_versor(self.particle.velocity) * lowest_time_to_collision)
        else:
            lowest_dist_to_collision, lowest_collision_segment, next_pos = \
                self.dist_intersection(remaining_dist, intersection_points, traveled_path)
            lowest_time_to_collision = lowest_dist_to_collision / norm(self.particle.velocity)

        return lowest_time_to_collision, lowest_dist_to_collision, lowest_collision_segment, next_pos


    def time_intersection(self, remaining_time, intersection_points, traveled_path):
        lowest_time_to_collision = remaining_time
        next_pos = None
        lowest_collision_segment = None
        for intersection_point, collision_element in intersection_points:
            time_to_collision = self._time_to_collision(intersection_point, traveled_path)
            if time_to_collision < lowest_time_to_collision:
                lowest_time_to_collision = time_to_collision
                lowest_collision_segment = collision_element
                next_pos = point_to_vec(intersection_point)
        lowest_time_to_collision *= TIME_PRECISION
        if not next_pos:
            next_pos = self.particle.calc_next_position(lowest_time_to_collision, None, self.check_condition)
        return lowest_time_to_collision, lowest_collision_segment, next_pos


    def dist_intersection(self, remaining_dist, intersection_points, traveled_path):
        lowest_dist_to_collision = remaining_dist
        next_pos = None
        lowest_collision_segment = None
        for intersection_point, collision_element in intersection_points:
            dist_to_collision = norm(traveled_path[0] - intersection_point)
            if dist_to_collision < lowest_dist_to_collision:
                lowest_dist_to_collision = dist_to_collision
                lowest_collision_segment = collision_element
                next_pos = point_to_vec(intersection_point)
        lowest_dist_to_collision *= DIST_PRECISION
        if not next_pos:
            next_pos = self.particle.calc_next_position(None, lowest_dist_to_collision, self.check_condition)
        return lowest_dist_to_collision, lowest_collision_segment, next_pos


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


    def _calc_stop_conditions(self, remaining_time, remaining_dist, count_loop):
        """
        Calculate stop conditions
        :param remaining_time: defined remaining_time
        :param remaining_dist: defined remaining_dist
        :param count_loop: check if loop is still executing
        :return: bool
        """
        if self.check_condition == "time":
            stop_conditions = np.isclose(
                remaining_time, 0, atol=10 ** (-self.significant_digits_time)) or remaining_time < 0
        else:
            stop_conditions = np.isclose(
                remaining_dist, 0, atol=10 ** (-self.significant_digits_dist)) or remaining_dist < 0
        if count_loop >= MAX_LOOP:
            print("Max loop iterations")
            return stop_conditions, True
        return stop_conditions, False


    def _time_to_collision(self, position: Point2, path: Segment2) -> float:
        """
        Calculate time until collision happen. Consider uniform particle movement in delta_t

        :param position: particle possible next position
        :param path: estimated particle path
        :return: time until collision point
        """
        pos = Segment2(path[0], position)
        return np.sqrt(float(pos.squared_length() / self.particle.velocity.squared_length()))


    def save_particle_data(self):
        """
        Save specific particle positions (if TEST is equal True)

        :return: None
        """
        if TEST and self.particle.id == followed_particle_id:
            particle_pos = Point2(self.particle.position.x(), self.particle.position.y())
            self.particle.positions.append(particle_pos)
