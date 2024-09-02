import itertools
from typing import Union
from skgeom.draw import draw
from skgeom import Bbox2, Vector2, Segment2
from utils.probabilistic_operations import random_vec
from scipy.constants import elementary_charge, electron_mass
from utils.complementary_operations import mirror, calc_versor


STOP_CONDITION = 10


class Particle:
    id_iter = itertools.count()


    def __init__(self, density: float, drift_method: str, effective_mass: float, fermi_velocity: float, position=None):
        """
        Particle that will be simulated

        :param density: density of particles [electrons/particle]
        :param drift_method: defines how drift velocity is calculated
        :param effective_mass: effective electron mass according material specifications [kg]
        :param fermi_velocity: scalar Fermi velocity [m/s]
        :param position: initial position (if desired to define nonrandom initial position)
        """
        self.id = next(Particle.id_iter)
        self.density = density
        self.charge = (-1) * elementary_charge * self.density
        self.effective_mass = effective_mass
        self.mass = self.effective_mass * self.density * electron_mass
        self.acceleration = None
        self.scalar_fermi_velocity = fermi_velocity
        self.fermi_velocity = None
        self.velocity = None
        self.position = position
        self.positions = list()
        self.drift_method = drift_method


    def set_init_position(self, bbox: Bbox2):
        """
        Set particle initial position into box bbox

        :param bbox: box size responsible to limit initial possible positions
        :return: None
        """
        min_range = (bbox.xmin(), bbox.ymin())
        max_range = (bbox.xmax(), bbox.ymax())
        self.position = Vector2(*random_vec(min_value=min_range, max_value=max_range, is_normalized=False))


    def set_velocity(self):
        """
        Set particle random Fermi velocity

        :return: None
        """
        self.fermi_velocity = Vector2(*random_vec()) * self.scalar_fermi_velocity
        self.velocity = self.fermi_velocity


    def calc_acceleration(self, electric_field: Vector2) -> Vector2:
        """
        Calculate particle acceleration

        :param electric_field: applied electric field in particle position
        :return: possible new acceleration
        """
        acceleration = electric_field * (self.charge / self.mass)
        return acceleration


    def calc_next_position(self, delta_t: Union[float, None], delta_s: Union[float, None], check_condition) -> Vector2:
        """
        Calculate next particle position according kinematic equations for uniformly varied motion

        :param delta_t: time interval
        :param delta_s: distance interval
        :param check_condition: id to identify which method should be used to calculate next position
        :return: line segment that connect initial and final particle positions
        """
        if check_condition == "time":
            next_pos = self.position + self.velocity * delta_t
        else:
            next_pos = self.position + calc_versor(self.velocity) * delta_s

        return next_pos


    def mirror_particle(self, normal_vec):
        """
        Mirror particle when collide with segment

        :param normal_vec: collided segment normal vector
        :return: mirrored velocity
        """
        self.velocity = mirror(self.velocity, normal_vec)


    def calc_drift_velocity(self, relax_time, electric_field: Vector2, mobility: float) -> None:
        """
        Calculate particle drift velocity

        :param relax_time: material relaxation time
        :param electric_field: applied electric field
        :param mobility: electronic material mobility
        """
        if self.drift_method == 'relax':
            self.velocity += self.charge * relax_time * electric_field / self.mass
        else:
            self.velocity += -mobility * electric_field


    def plot_traveled_path(self):
        pos0 = self.positions[0]
        for pos in self.positions[1:]:
            draw(Segment2(pos0, pos))
            pos0 = pos
