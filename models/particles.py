from skgeom import Vector2

from models.particle import Particle
from models.electric_field import ElectricField
from utils.probabilistic_operations import random_vec
from utils.complementary_operations import vec_to_point
from scipy.constants import elementary_charge, electron_mass


class Particles:
    def __init__(
            self,
            effective_mass: float,
            fermi_velocity: float,
            n_particles: int,
            density: float
    ):
        self.density = density
        self.effective_mass = effective_mass
        self.scalar_fermi_velocity = fermi_velocity
        self.charge = (-1) * elementary_charge * self.density
        self.particles_list = self.create_particles(n_particles)
        self.mass = self.effective_mass * self.density * electron_mass


    def set_velocity_fermi(self, particle: Particle = None):
        """
        Set particle random Fermi velocity

        :return: None
        """
        if not particle:
            for particle in self.particles_list:
                particle.velocity_fermi = Vector2(*random_vec()) * self.scalar_fermi_velocity
        else:
            particle.velocity_fermi = Vector2(*random_vec()) * self.scalar_fermi_velocity


    def set_velocity_drift(self, relax_time, electric_field: ElectricField, particle: Particle):
        """
        Calculate particle drift velocity

        :param relax_time: material relaxation time
        :param electric_field: applied electric field
        :param particle:
        :return: drift velocity
        """
        particle.velocity_drift = self.charge * relax_time * electric_field.vector(particle.position) / self.mass


    def set_velocity_total(self, particle: Particle = None):
        """
        Set particle total velocity as a composition of Fermi velocity and drift velocity
        :return:
        """
        if not particle:
            for particle in self.particles_list:
                particle.velocity_total = particle.velocity_fermi + particle.velocity_drift
        else:
            particle.velocity_total = particle.velocity_fermi + particle.velocity_drift


    def set_acceleration(self, electric_field: Vector2, particle: Particle):
        """
        Calculate particle acceleration

        :param electric_field: applied electric field in particle position
        :param particle:
        :return: possible new acceleration
        """
        particle.acceleration = electric_field * (self.charge / self.mass)


    def get_all_pos(self):
        """
        Get all particle positions
        :return: list of positions
        """
        return [vec_to_point(particle.position) for particle in self.particles_list]


    @staticmethod
    def create_particles(number_of_particles) -> list:
        """
        Create particles according specified number or according to CPU cores

        :param number_of_particles: number of particles to be created (if specified)
        :return: list of particles
        """
        particles = list()
        for _ in range(number_of_particles):
            particles.append(Particle())
        return particles
