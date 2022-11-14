from skgeom.draw import draw
from material import Material
from utils.probabilistic_operations import random_vec
from scipy.constants import electron_mass, elementary_charge
from skgeom import Bbox2, Vector2, Point2, Segment2
from utils.complementary_operations import mirror, vec_to_point


class Particle:
    def __init__(self, density: float, effective_mass: float, size: float, fermi_velocity: float, position=None):
        self.charge = elementary_charge * density
        self.mass = electron_mass * effective_mass * density
        self.size = size
        self.acceleration = None
        self.scalar_fermi_velocity = fermi_velocity
        self.fermi_velocity = random_vec() * self.scalar_fermi_velocity
        self.velocity = self.fermi_velocity
        self.position = position
        self.positions = list()
        self.travelling_time = 0


    def set_init_position(self, bbox: Bbox2):
        """
        Set particle initial position into box bbox

        :param bbox: box size
        :return: None
        """
        min_range = (bbox.xmin(), bbox.ymin())
        max_range = (bbox.xmax(), bbox.ymax())
        self.position = random_vec(min_value=min_range, max_value=max_range, is_normalized=False)


    def set_fermi_velocity(self):
        """
        Set particle random Fermi velocity

        :return: None
        """
        self.fermi_velocity = random_vec() * self.scalar_fermi_velocity


    def calc_acceleration(self, electric_field: Vector2):
        """
        Calculate particle acceleration

        :param electric_field: applied electric field in particle position
        :return: None
        """
        self.acceleration = electric_field * (-1 * self.charge / self.mass)


    def calc_velocity(self, electric_field: Vector2, delta_t: float):
        """
        Calculate particle velocity after delta_t

        :param electric_field: applied electric field in particle position
        :param delta_t: time interval
        :return: None
        """
        self.calc_acceleration(electric_field)
        self.velocity += self.acceleration * delta_t


    def calc_position(self, delta_t: float):
        """
        Calculate next particle position according cinematic equations for uniformly varied motion

        :param delta_t: time interval
        :return: line segment that connect initial and final particle positions
        """
        next_pos = self.position + self.velocity * delta_t + (self.acceleration * delta_t ** 2) / 2
        p_0 = vec_to_point(self.position)
        p_1 = vec_to_point(next_pos)
        path = Segment2(p_0, p_1)
        return path


    def time_to_collision(self, position: Point2, path: Segment2):
        """
        Calculate time interval until defined collision

        :param position: particle position
        :param path: crossed line segment
        :return: time until collision
        """
        pos = Segment2(path[0], position)
        return float(pos.squared_length() / self.velocity.squared_length()) ** (1 / 2)


    def move(self, normal_vec: Vector2, delta_t: float, electric_field: Vector2, relaxation: bool):
        """
        Move particle according time interval

        :param normal_vec: boundary normal vector
        :param delta_t: time interval
        :param electric_field: particle applied electric field
        :param relaxation: relaxation event
        :return: None
        """
        self.position = self.position + self.velocity * delta_t
        self.velocity = mirror(self.velocity, normal_vec)
        if relaxation:
            self.set_fermi_velocity()
            self.calc_velocity(electric_field, delta_t)


if __name__ == '__main__':
    particle = Particle(1, 1, 1, 10)
    mat = Material(1, 1, 1)
    e_field = Vector2(1, 0)
    particle.calc_velocity(e_field, 1)
    box = Bbox2(1, 2, 3, 4)
    particle.set_init_position(box)
    print('ei')
