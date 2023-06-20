import matplotlib
import numpy as np

from models.topology import Topology
from models.particle import Particle
from models.material import Material
from models.particles import Particles
from skgeom import Vector2, Point2, Segment2
from models.electric_field import ElectricField
from utils.probabilistic_operations import random_vec
from utils.complementary_operations import vec_to_point, point_to_vec

TEST = False
TIME_PRECISION = 0.99
matplotlib.use('TkAgg')


class System:
    def __init__(self, data: dict):
        """
        Create system to be simulated (topology + particles + materials + etc.)
        :param data: number of particles (defined by the number of cores if not defined)
        """

        self.material = Material(**data['material'])
        self.particles = Particles(
            self.material.effective_mass,
            self.material.scalar_fermi_velocity,
            **data['particles']
        )
        self.topology = self.chose_topology(data['geometry'])
        self.electric_field = ElectricField()

        self.total_macro_particles = len(self.particles.particles_list) * self.particles.density

        self.particles_counter = 0
        self.simulated_time = 0
        self.time_steps_count = 0
        self.collisions_count = 0


    def reset_performance_params(self):
        """
        Reset all performance parameters for the next simulation
        :return:
        """
        self.particles_counter = 0
        self.simulated_time = 0
        self.time_steps_count = 0
        self.collisions_count = 0


    def set_particles_parameters(self):
        """
        Set particles initial position and Fermi velocity

        :return: None
        """
        def set_pos(min_r, max_r):
            pos = Vector2(*random_vec(min_value=min_r, max_value=max_r, is_normalized=False))
            i_pos = vec_to_point(pos)
            _, l_dist = self.topology.get_closer_segment(i_pos)
            return pos, i_pos, l_dist

        min_range = (self.topology.bbox.xmin(), self.topology.bbox.ymin())
        max_range = (self.topology.bbox.xmax(), self.topology.bbox.ymax())

        for particle in self.particles.particles_list:
            particle.position, init_pos, lowest_dist = set_pos(min_range, max_range)

            while (not self.topology.contains(init_pos)) or (lowest_dist < self.topology.scale / 20):
                particle.position, init_pos, lowest_dist = set_pos(min_range, max_range)


    def check_current_segment_collision(self, particle: Particle, closest_collision_segment, segment_normal_vec):
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
            self.particles_counter += self.particles.density
            return True, 'reverse'
        elif collided_element and collided_element in self.topology.current_computing_elements['reverse']:
            self.particles_counter -= self.particles.density
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
        particle.position = point_to_vec(pos)


    def calc_closer_intersection(self, remaining_time: float, particle: Particle) -> (float, Segment2, Vector2):
        """
        Define closer intersection between particle path and geometries boundaries

        :param remaining_time: time until scattering process
        :param particle: simulated particle
        :return lowest_time_to_collision: the lowest time to collision
        :return lowest_collision_segment: collided segment
        :return next_pos: particle next position
        """
        lowest_time_to_collision = remaining_time
        next_pos = None
        lowest_collision_segment = None
        intersection_points = self.topology.intersection_points(particle.travelled_path)
        for intersection_point, collision_element in intersection_points:
            time_to_collision = self._time_to_collision(particle, intersection_point)
            if time_to_collision < lowest_time_to_collision:
                lowest_time_to_collision = time_to_collision * TIME_PRECISION
                lowest_collision_segment = collision_element
                next_pos = intersection_point

        if not next_pos:
            next_pos = particle.position + particle.velocity_total * lowest_time_to_collision

        return lowest_time_to_collision, lowest_collision_segment, next_pos


    def set_particle_velocity_drift(self, particle):
        """
        Set particle drift velocity
        :param particle: particle to be changed
        :return:
        """
        self.particles.set_velocity_drift(self.material.relax_time, self.electric_field, particle)


    def scatter_particles(self, particle: Particle = None):
        """
        Scatter all particles or specific particle
        :param particle: specified particle
        :return:
        """
        self.particles.set_velocity_fermi(particle)
        self.particles.set_velocity_total(particle)


    @staticmethod
    def _time_to_collision(particle: Particle, position: Point2) -> float:
        """
        Calculate time until collision happen. Consider uniform particle movement in delta_t

        :param position: particle possible next position
        :return: time until collision point
        """
        pos = Segment2(particle.travelled_path[0], position)
        return np.sqrt(float(pos.squared_length() / particle.velocity_total.squared_length()))


    @staticmethod
    def chose_topology(geometry_dict) -> Topology:
        """
        Chose how to set topology
        :param geometry_dict: input dictionary
        :return: Topology instance
        """
        geometry_type = geometry_dict['input_style']
        del geometry_dict['input_style']
        if geometry_type == 'file':
            return Topology.from_file(**geometry_dict)
        else:
            return Topology.from_points(**geometry_dict)
