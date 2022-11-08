import random
from skgeom import Vector2
from utils.complementary_operations import norm


def random_vec(shape=2, min_value: tuple = (-1, -1), max_value: tuple = (1, 1), is_normalized: bool = True) -> Vector2:
    random_numbers = list()
    for pos in range(shape):
        random_numbers.append(random.uniform(min_value[pos], max_value[pos]))
    rand_vec = Vector2(*random_numbers)
    if is_normalized:
        vec_norm = norm(rand_vec)
        rand_vec = rand_vec / vec_norm
    return rand_vec
