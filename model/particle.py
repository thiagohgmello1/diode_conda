from utils.random_gen import random_vec
from utils.complementary_operations import mirror, vec_to_point
from skgeom import Bbox2, Vector2, Point2, Segment2
from skgeom.draw import draw
from material import Material


class Particle:
    def __init__(self, charge: float, mass: float, size: float, fermi_velocity: float, position=None):
        self.charge = charge
        self.mass = mass
        self.size = size
        self.scalar_fermi_velocity = fermi_velocity
        self.fermi_velocity = random_vec() * self.scalar_fermi_velocity
        self.velocity = None
        self.position = position
        self.positions = list()


    def set_fermi_velocity(self):
        self.fermi_velocity = random_vec() * self.scalar_fermi_velocity


    def calc_velocity(self, electric_field: Vector2, material: Material):
        avg_vel = electric_field * (-1 * self.charge * material.relax_time / (self.mass * material.mass_scale))
        self.velocity = self.fermi_velocity + avg_vel


    def set_init_position(self, bbox: Bbox2):
        min_range = (bbox.xmin(), bbox.ymin())
        max_range = (bbox.xmax(), bbox.ymax())
        self.position = random_vec(min_value=min_range, max_value=max_range, is_normalized=False)


    def calc_next_pos(self, time: float):
        next_pos = self.position + self.velocity * time
        p_0 = vec_to_point(self.position)
        p_1 = vec_to_point(next_pos)
        path = Segment2(p_0, p_1)
        return path


    def time_to_collision(self, position: Point2, path: Segment2):
        pos = Segment2(path[0], position)
        return float(pos.squared_length() / self.velocity.squared_length()) ** (1 / 2)


    def move(self, normal_vec: Vector2, time: float, electric_field: Vector2, material: Material):
        self.position = self.position + self.velocity * time
        self.velocity = mirror(self.velocity, normal_vec)
        # self.calc_velocity(electric_field, material)


if __name__ == '__main__':
    particle = Particle(1, 1, 1, 10)
    mat = Material(1, 1)
    e_field = Vector2(1, 0)
    particle.calc_velocity(e_field, mat)
    box = Bbox2(1, 2, 3, 4)
    particle.set_init_position(box)
    print('ei')
