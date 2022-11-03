import random

import numpy as np


def random_vec(shape=2, min_value: tuple = (-1, -1), max_value: tuple = (1, 1), is_normalized: bool = True):
    random_numbers = list()
    for pos in range(shape):
        random_numbers.append(random.uniform(min_value[pos], max_value[pos]))
    rand_vec = np.array(random_numbers)
    if is_normalized:
        norm = np.linalg.norm(rand_vec)
        rand_vec /= norm
    return rand_vec
