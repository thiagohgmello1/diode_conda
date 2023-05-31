from utils.probabilistic_operations import random_vec
from scipy.constants import elementary_charge, electron_mass
from skgeom import Bbox2, Vector2, Point2, Segment2
from utils.complementary_operations import mirror, vec_to_point


STOP_CONDITION = 10


class Particle:
    def __init__(self, density: float, effective_mass: float, fermi_velocity: float, position=None):
        """
        Particle that will be simulated

        :param density: density of particles [electrons/particle]
        :param effective_mass: effective electron mass according material specifications [kg]
        :param fermi_velocity: scalar Fermi velocity [m/s]
        :param position: initial position (if desired to define nonrandom initial position)
        """
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


    def calc_next_position(self, delta_t: float) -> Segment2:
        """
        Calculate next particle position according cinematic equations for uniformly varied motion

        :param delta_t: time interval
        :return: next possible position and line segment that connect initial and final particle positions
        """
        next_pos = self.position + self.velocity * delta_t
        p_0 = vec_to_point(self.position)
        p_1 = vec_to_point(next_pos)
        path = Segment2(p_0, p_1)
        return path


    def mirror_particle(self, normal_vec):
        """
        Mirror particle when collide with segment

        :param normal_vec: collided segment normal vector
        :return: mirrored velocity
        """
        self.velocity = mirror(self.velocity, normal_vec)


    @staticmethod
    def calc_drift_velocity(mobility: float, electric_field: Vector2) -> Vector2:
        """
        Calculate particle drift velocity

        :param mobility: material mobility
        :param electric_field: applied electric field
        :return: drift velocity
        """
        return mobility * electric_field
