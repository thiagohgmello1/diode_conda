from itertools import count
from dataclasses import dataclass, field
from skgeom import Vector2, Segment2
from utils.complementary_operations import mirror, vec_to_point


@dataclass
class Particle:
    id: int = field(default_factory=count().__next__, init=False, repr=True)
    acceleration: Vector2 = field(default=Vector2(0, 0), init=False)
    velocity_total: Vector2 = field(default=Vector2(0, 0), init=False)
    velocity_fermi: Vector2 = field(default=Vector2(0, 0), init=False)
    velocity_drift: Vector2 = field(default=Vector2(0, 0), init=False)
    position: Vector2 = field(default=None, init=False)
    travelled_path: Segment2 = field(default=None, init=False)


    def mirror_particle(self, normal_vec):
        """
        Mirror particle when collide with segment

        :param normal_vec: collided segment normal vector
        :return: mirrored velocity
        """
        self.velocity_total = mirror(self.velocity_total, normal_vec)


    def calc_next_pos_uniform(self, delta_t: float):
        """
        Calculate next particle position according kinematic equations for uniform motion

        :param delta_t: time interval
        """
        next_pos = self.position + self.velocity_total * delta_t
        p_0 = vec_to_point(self.position)
        p_1 = vec_to_point(next_pos)
        path = Segment2(p_0, p_1)
        self.travelled_path = path
