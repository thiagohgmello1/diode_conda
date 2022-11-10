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


    def calc_acceleration(self, electric_field: Vector2):
        self.acceleration = electric_field * (-1 * self.charge / self.mass)


    def set_fermi_velocity(self):
        self.fermi_velocity = random_vec() * self.scalar_fermi_velocity


    def calc_velocity(self, electric_field: Vector2, delta_t: float):
        self.calc_acceleration(electric_field)
        self.velocity += self.acceleration * delta_t


    def set_init_position(self, bbox: Bbox2):
        min_range = (bbox.xmin(), bbox.ymin())
        max_range = (bbox.xmax(), bbox.ymax())
        self.position = random_vec(min_value=min_range, max_value=max_range, is_normalized=False)


    def calc_position(self, time: float):
        next_pos = self.position + self.velocity * time
        p_0 = vec_to_point(self.position)
        p_1 = vec_to_point(next_pos)
        path = Segment2(p_0, p_1)
        return path


    def time_to_collision(self, position: Point2, path: Segment2):
        pos = Segment2(path[0], position)
        return float(pos.squared_length() / self.velocity.squared_length()) ** (1 / 2)


    def move(self, normal_vec: Vector2, delta_t: float, electric_field: Vector2, relaxation: bool):
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