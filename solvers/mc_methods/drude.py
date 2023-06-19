import numpy as np

from solvers.mc_methods.method import Method
from scipy.constants import elementary_charge
from utils.complementary_operations import vec_to_point, calc_normal


class Drude(Method):
    def __init__(self, system):
        super().__init__(system)


    def simulate(self, particle):
        self.system.set_particle_velocity_drift(particle)
        particle.velocity_total += particle.velocity_drift
        remaining_time = self.system.material.relax_time

        stop_conditions = np.isclose(
            remaining_time, 0, atol=10 ** (-self.significant_digits_time)) or remaining_time < 0
        while not stop_conditions:
            particle.calc_next_pos_uniform(remaining_time)
            lowest_time_to_collision, closest_collision_segment, next_pos = self.system.calc_closer_intersection(
                remaining_time, particle
            )
            particle_p0 = vec_to_point(particle.position)
            particle.position = next_pos

            if not self.system.topology.contains(vec_to_point(particle.position)):
                raise Exception('Particle is outside geometry')

            collision_normal_vec = calc_normal(closest_collision_segment, particle_p0)
            self.system.check_current_segment_collision(particle, closest_collision_segment, collision_normal_vec)
            remaining_time -= lowest_time_to_collision
            stop_conditions = np.isclose(
                remaining_time, 0, atol=10 ** (-self.significant_digits_time)) or remaining_time < 0


    def calc_current(self):
        """
        Calculate total current

        :return: calculated current
        """
        carrier_c = self.system.material.carrier_concentration
        current = (
                          carrier_c * self.system.topology.area * elementary_charge * self.system.particles_counter
                  ) / (self.system.simulated_time * self.system.total_macro_particles)
        return current
