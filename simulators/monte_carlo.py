import time
import matplotlib
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

from skgeom.draw import draw
from datetime import datetime
from math import log10, floor
from model.particle import Particle
from model.topology import Topology
from model.material import Material
from skgeom import Vector2, Point2, Segment2
from utils.comparable_methods import drude_analytical_model
from scipy.constants import c, elementary_charge, electron_mass
from utils.probabilistic_operations import decision, random_number
from utils.complementary_operations import vec_to_point, calc_normal, create_segments

TEST = False
SIGNIFICANT_DIGITS = 6
RESOLUTION = 1
BREAK_MAX = 1000
SCALE = 10
matplotlib.use('TkAgg')


class System:
    def __init__(
            self,
            particle: Particle,
            topology: Topology,
            material: Material,
            electric_field: Vector2,
            max_simulations: float = np.inf,
            number_of_particles: int = None,
            max_time_simulation: float = np.inf,
            max_collisions: float = np.inf,
            max_time_steps: float = np.inf
    ):
        """
        Create system to be simulated (topology + particles + materials + etc.)

        :param particle: particle model to be simulated
        :param topology: desired topology
        :param material: material to topology
        :param electric_field: defined or calculated electric field created by applied voltage in topology [V/m]
        :param max_simulations: number of simulated particles
        :param number_of_particles: number of particles (defined by the number of cores if not defined)
        :param max_time_simulation: maximum time to simulation [s]
        :param max_collisions: defined maximum accepted collisions. Stop criteria
        :param max_time_steps: defined maximum time steps. Stop criteria [s]
        """
        self.max_simulations = max_simulations
        self.particles = self.create_particles(particle, number_of_particles)
        self.topology = topology
        self.material = material
        self.e_field = electric_field
        self.relax_time = self.material.relax_time

        self.total_macro_particles = 0
        self.particle_counter = 0

        self.max_time_simulation = max_time_simulation
        self.simulated_time = list()

        self.max_collisions = max_collisions
        self.collisions = 0

        self.max_time_steps = max_time_steps
        self.time_steps = 0

        self.simulations_counter = 0
        self.significant_digits = -(int(floor(log10(self.material.relax_time)))) + SIGNIFICANT_DIGITS - 1


    @staticmethod
    def create_particles(particle_model, number_of_particles) -> list:
        """
        Create particles according specified number or according CPU cores

        :param particle_model: model of desired particle
        :param number_of_particles: number of particles to be created (if specified)
        :return: list of particles
        """
        particles = list()
        if not number_of_particles:
            number_of_particles = multiprocessing.cpu_count()
        density = particle_model.density
        effective_mass = particle_model.effective_mass
        fermi_velocity = particle_model.scalar_fermi_velocity
        for _ in range(number_of_particles):
            particles.append(Particle(density, effective_mass, fermi_velocity))
        return particles


    def _set_time_step(self, time_step) -> float:
        """
        Set simulation time step

        :param time_step: input time step representing a fraction of relaxation time (if desired)
        :return: simulation time step
        """
        if not time_step:
            max_e_field = np.sqrt(float(self.e_field.squared_length()))
            min_dimension = min(
                [
                    self.topology.bbox.xmax() - self.topology.bbox.xmin(),
                    self.topology.bbox.ymax() - self.topology.bbox.ymin()
                ]
            )
            fermi_min = min_dimension / self.material.scalar_fermi_velocity
            drift_min = np.sqrt(min_dimension * self.particles[0].mass / (max_e_field * abs(self.particles[0].charge)))
            return float(1 / SCALE * min(fermi_min, drift_min, self.relax_time))
        return self.relax_time / time_step


    def set_particle_parameters(self, particle):
        """
        Set particle initial parameters

        :param particle: particle to be started
        :return: None
        """
        particle.set_velocity()
        particle.set_init_position(self.topology.bbox)
        init_pos = vec_to_point(particle.position)

        while not self.topology.contains(init_pos):
            particle.set_init_position(self.topology.bbox)
            init_pos = vec_to_point(particle.position)
        if TEST:
            particle.positions.append([init_pos])


    def set_particles_parameters(self):
        """
        Set particles initial parameters

        :return: None
        """
        for particle in self.particles:
            self.set_particle_parameters(particle)


    def _relaxation_event(
            self,
            delta_t: float,
            cumulative_time: float,
            collision_segment: Segment2
    ) -> (str, float, Segment2):
        """
        Define if relaxation event occurs in time interval

        :param delta_t: time interval
        :param cumulative_time: cumulative time since the last collision
        :param collision_segment: segment where particle will collide
        :return relaxation: string indicating if there is or there is not relaxation event
        :return delta_t: time interval
        :return collision_segment: segment where particle will collide if there is no relaxation event
        """
        if cumulative_time >= self.relax_time:
            return 'relax', 0, None
        else:
            relax_probability = 1 - np.exp(-cumulative_time / self.relax_time)
        relaxation = decision(relax_probability)

        if relaxation:
            delta_t = random_number(0, delta_t)
            collision_segment = None
            relaxation = "collide"
        else:
            relaxation = "not collide"

        return relaxation, delta_t, collision_segment


    def simulate(self, model):
        """
        Simulate complete system. Simulation can be executed in GPU or single-multi core CPU

        :param model: desired model to simulate. Must be a method callback (ex.: simulate_drude method)
        :return: None
        """
        while self._calc_stop_conditions():
            self.total_macro_particles += self.particles[0].density
            self.set_particle_parameters(self.particles[0])
            self.simulations_counter += 1
            model(self.particles[0])
            progress_bar(self.collisions, self.max_collisions)
        print('\n')


    def _calc_stop_conditions(self) -> bool:
        """
        Calculate if any stop condition was met

        :return: stop condition
        """
        time_steps_condition = self.time_steps < self.max_time_steps
        collisions_condition = self.collisions < self.max_collisions
        simulations_condition = self.simulations_counter < self.max_simulations

        return time_steps_condition and collisions_condition and simulations_condition


    def simulate_drude(self, particle: Particle):
        """
        Simulate Drude event

        :param particle: particle to be simulated
        :return: None
        """
        particle.acceleration = particle.calc_acceleration(self.e_field)
        particle.velocity += particle.acceleration * self.relax_time

        remaining_time = self.relax_time
        stop_conditions = np.isclose(remaining_time, 0, atol=0)
        break_counter = 0

        while (not stop_conditions) and (break_counter < BREAK_MAX):
            break_counter += 1
            if remaining_time < 0:
                print(f'ERROR: {remaining_time}')
            traveled_path = particle.calc_next_position(remaining_time)
            intersection_points = self.topology.intersection_points(traveled_path)
            lowest_time_to_collision, closest_collision_segment = self._calc_closer_intersection(
                remaining_time, particle.velocity, intersection_points, traveled_path
            )
            particle.position = particle.position + particle.velocity * lowest_time_to_collision
            particle_position = Point2(particle.position.x(), particle.position.y())
            segment_normal_vec = calc_normal(closest_collision_segment, particle_position)
            particle.mirror_particle(segment_normal_vec)
            remaining_time -= lowest_time_to_collision
            remaining_time = round(remaining_time, self.significant_digits)

            if segment_normal_vec:
                self.collisions += 1
                current_collision = self.particle_computation(closest_collision_segment, particle.density)
                if current_collision:
                    stop_conditions = True
                    self.save_particle_data(particle)
                    continue
            self.time_steps += 1
            self.save_particle_data(particle)
            stop_conditions = np.isclose(remaining_time, 0, atol=0, rtol=self.significant_digits)
        if break_counter >= BREAK_MAX:
            print(f'ERROR: E = {self.e_field}')


    def particle_computation(self, collided_element: Segment2, particle_density: float) -> bool:
        """
        Compute collisions in current elements

        :param collided_element: collided geometry segment
        :param particle_density: particle density
        :return: boolean indicating if there was rectification
        """
        if collided_element and collided_element in self.topology.current_computing_elements['direct']:
            self.particle_counter += particle_density
            return True
        elif collided_element and collided_element in self.topology.current_computing_elements['reverse']:
            self.particle_counter -= particle_density
            return True
        return False


    def _calc_closer_intersection(
            self,
            remaining_time,
            particle_velocity: Vector2,
            intersection_points: list[Point2],
            traveled_path: Segment2
    ) -> (float, Segment2):
        """
        Define closer intersection between particle path and geometries boundaries

        :param intersection_points: list of points where collision can occur
        :param particle_velocity: possible particle velocity
        :param traveled_path: corresponding particle path
        :return lowest_time_to_collision: lowest time to collision
        :return lowest_collision_segment: collided segment
        """
        lowest_time_to_collision = remaining_time
        lowest_collision_segment = None
        for intersection_point, collision_element in intersection_points:
            time_to_collision = self._time_to_collision(particle_velocity, intersection_point, traveled_path)
            if time_to_collision < lowest_time_to_collision:
                lowest_time_to_collision = time_to_collision * RESOLUTION
                lowest_collision_segment = collision_element

        return lowest_time_to_collision, lowest_collision_segment


    def cal_current(self):
        """
        Calculate total current

        :return: calculated currents
        """
        carrier_concentration = self.material.carrier_concentration
        current = (carrier_concentration * self.topology.area * elementary_charge * self.particle_counter) /\
                  (self.total_macro_particles * self.max_time_simulation)
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


    def save_particle_data(self, particle):
        """
        Save particle positions (if TEST is equal True)

        :param particle: desired particle to save data
        :return: None
        """
        if TEST:
            particle_pos = Point2(particle.position.x(), particle.position.y())
            particle.positions[self.simulations_counter - 1].append(particle_pos)


def draw_behaviour(desired_sys, simulations):
    """
    Draw some achieved particle positions

    :param desired_sys: system to be drawn
    :return: None
    """
    if TEST:
        draw(desired_sys.topology.topologies)
        segments_to_print = list()
        for simulation in simulations:
            segments_to_print.append(create_segments(desired_sys.particles[0].positions[simulation]))
        draw(segments_to_print)
        plt.savefig('../outputs/diode.png')


def save_current(current_file: str, meas_current: float, simulated_geometry: str, electric_field: list):
    """
    Save calculated current

    :param current_file: output file
    :param meas_current: calculated current
    :param simulated_geometry: .svg simulated file
    :return: None
    """
    now = datetime.now()
    date_string = now.strftime("%d/%m/%Y %H:%M:%S")
    with open(current_file, 'a') as f:
        string_to_be_saved = f'{date_string};{simulated_geometry};{meas_current};{electric_field}\n'
        f.write(string_to_be_saved)


def progress_bar(progress, total):
    percent = 100 * (progress / total)
    bar = '█' * int(percent) + '-' * (100 - int(percent))
    print(f'\r|{bar}| {percent:.2f}%', end='\r')


if __name__ == '__main__':
    f_velocity = c / 300
    MFPL = 500e-9
    carrier_c = 7.2e15
    thickness = 300e-9
    gate_voltage = 10
    geometry = '../tests/diode7.svg'
    mat = Material(
        mean_free_path=MFPL,
        scalar_fermi_velocity=f_velocity,
        permittivity=3.9,
        substrate_thickness=thickness,
        gate_voltage=gate_voltage
    )
    particle_m = Particle(density=120, effective_mass=mat.effective_mass, fermi_velocity=mat.scalar_fermi_velocity)
    pol = Topology.from_file(geometry, 1e-7)
    e_field = Vector2(-0.5, 0) / (pol.bbox.xmax() - pol.bbox.xmin())
    system = System(
        particle=particle_m,
        topology=pol,
        material=mat,
        electric_field=e_field,
        max_collisions=100000,
        max_time_simulation=mat.relax_time
    )
    exec_time = time.time()
    system.simulate(system.simulate_drude)
    exec_time = time.time() - exec_time

    drude_analytical_current = drude_analytical_model(
        width=float(pol.bbox.ymax() - pol.bbox.ymin()),
        relax_time=mat.relax_time,
        carrier_concentration=mat.carrier_concentration,
        effective_mass=mat.effective_mass * electron_mass,
        e_field=e_field
    )
    simulation_current = system.cal_current()
    save_current('../outputs/currents.csv', simulation_current, geometry, e_field)
    print(f'Time steps: {system.time_steps}')
    print(f'Collisions: {system.collisions}')
    print(f'Current:{simulation_current}')
    print(f'Drude current: {drude_analytical_current}')
    print(f'Execution time: {exec_time}')
    draw_behaviour(system, [i for i in range(1000, 1100)])
