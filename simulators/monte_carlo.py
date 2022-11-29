import numpy as np

from skgeom.draw import draw
from scipy.constants import c
from model.particle import Particle
from model.topology import Topology
from model.material import Material
from skgeom import Vector2, Point2, Segment2
from utils.probabilistic_operations import decision, random_number
from utils.complementary_operations import vec_to_point, calc_normal, create_segments


TEST = True
SCALE = 10


class System:
    def __init__(
            self,
            particles: list[Particle],
            topology: Topology,
            material: Material,
            electric_field: Vector2,
            max_time_simulation: float = np.inf,
            max_collisions: float = np.inf,
            max_time_steps: float = np.inf,
            time_step: float = None
    ):
        """
        Create system to be simulated (topology + particles + materials + etc.)

        :param particles: particles to be simulated
        :param topology: desired topology
        :param material: material to topology
        :param electric_field: defined or calculated electric field created by applied voltage in topology [V/m]
        :param max_time_simulation: maximum time to simulation [s]
        :param max_collisions: defined maximum accepted collisions. Stop criteria
        :param max_time_steps: defined maximum time steps. Stop criteria [s]
        :param time_step: fraction of relaxation time to define discretization
        """
        self.particles = particles
        self.topology = topology
        self.material = material
        self.e_field = electric_field
        self.relax_time = self.material.relax_time
        self.time_step = self._set_time_step(time_step)

        self.collisions = 0
        self.time_steps = 0

        self.max_time_simulation = max_time_simulation
        self.max_collisions = max_collisions
        self.max_time_steps = max_time_steps


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
            particle.positions.append(init_pos)


    def _relaxation_event(
            self,
            delta_t: float,
            cumulative_time: float,
            collision_segment: Segment2
    ) -> (bool, float, Segment2):
        """
        Define if relaxation event occur in time interval

        :param delta_t: time interval
        :param cumulative_time: cumulative time since the last collision
        :param collision_segment: segment where particle will collide
        :return relaxation: boolean indicating if there is or there is not relaxation event
        :return delta_t: time interval
        :return collision_segment: segment where particle will collide if there is no relaxation event
        """
        relax_probability = 1 - np.exp(-cumulative_time / self.relax_time)
        relaxation = decision(relax_probability)

        if relaxation:
            delta_t = random_number(0, delta_t)
            collision_segment = None

        return relaxation, delta_t, collision_segment


    def _calc_particle_parameters(self, particle) -> Segment2:
        """
        Calculate particle mechanical parameters

        :param particle: particle to be set-up
        :return: possible particle traveled path
        """
        particle.acceleration = particle.calc_acceleration(self.e_field)
        particle.velocity += particle.acceleration * self.time_step
        particle_traveled_path = particle.calc_next_position(particle.velocity, self.time_step)
        return particle_traveled_path


    def simulate_drude(self, particle: Particle):
        """
        Simulate Drude event. Simulation can be executed in GPU or single-multithreading CPU

        :param particle: particle to be simulated
        :return: None
        """
        simulated_time = 0
        condition = self._stop_condition(simulated_time)
        cumulative_time = 0

        while not condition:
            traveled_path = self._calc_particle_parameters(particle)
            intersection_points = self.topology.intersection_points(traveled_path)
            lowest_time_to_collision, lowest_collision_segment = self._calc_closer_intersection(
                particle.velocity, intersection_points, traveled_path
            )
            cumulative_time += lowest_time_to_collision
            relaxation, lowest_time_to_collision, lowest_collision_segment = self._relaxation_event(
                lowest_time_to_collision, cumulative_time, lowest_collision_segment
            )
            particle_pos = particle.position + particle.velocity * lowest_time_to_collision
            particle_pos = Point2(particle_pos.x(), particle_pos.y())
            segment_normal_vec = calc_normal(lowest_collision_segment, particle_pos)
            lowest_time_to_collision = particle.move(
                segment_normal_vec, lowest_time_to_collision, relaxation, self.topology.contains
            )

            if relaxation:
                cumulative_time = 0
            elif segment_normal_vec:
                cumulative_time = 0
                self.collisions += 1

            simulated_time += lowest_time_to_collision
            self.time_steps += 1
            condition = self._stop_condition(simulated_time)
            if TEST:
                particle_pos = Point2(particle.position.x(), particle.position.y())
                particle.positions.append(particle_pos)


    def _stop_condition(self, simulated_time) -> bool:
        """
        Calculate if any stop condition was met

        :param simulated_time: total simulated time for specific particle
        :return: stop condition
        """
        time_steps_condition = self.time_steps > self.max_time_steps
        collisions_condition = self.collisions > self.max_collisions
        time_condition = np.isclose(simulated_time, self.max_time_simulation)

        return time_steps_condition or collisions_condition or time_condition


    def _calc_closer_intersection(
            self,
            particle_velocity,
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
        lowest_time_to_collision = self.time_step
        lowest_collision_segment = None
        for intersection_point, collision_element in intersection_points:
            time_to_collision = self._time_to_collision(particle_velocity, intersection_point, traveled_path)
            if time_to_collision < lowest_time_to_collision:
                lowest_time_to_collision = time_to_collision
                lowest_collision_segment = collision_element

        return lowest_time_to_collision, lowest_collision_segment


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
        return float(pos.squared_length() / particle_velocity.squared_length()) ** (1 / 2)


if __name__ == '__main__':
    f_velocity = c / 300
    MFPL = 200e-9
    carrier_c = 1.1e16
    mat = Material(mean_free_path=MFPL, scalar_fermi_velocity=f_velocity, carrier_concentration=carrier_c)
    particles_list = [Particle(density=10, effective_mass=mat.effective_mass, fermi_velocity=mat.scalar_fermi_velocity)]
    pol = Topology.from_file('../tests/test.svg', 1e-9)
    e_field = Vector2(-1, 0) / (pol.bbox.xmax() - pol.bbox.xmin())
    system = System(particles_list, pol, mat, e_field, max_collisions=20, max_time_steps=50)
    system.set_particles_parameters()
    system.simulate_drude(particles_list[0])

    draw(system.topology.topologies)
    if TEST:
        segments = create_segments(particles_list[0].positions)
        draw(segments)
    print('ei')
