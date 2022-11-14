import random
from skgeom import Vector2
from utils.complementary_operations import norm


def random_vec(shape=2, min_value: tuple = (-1, -1), max_value: tuple = (1, 1), is_normalized: bool = True) -> Vector2:
    """
    Calculate random vector

    :param shape: returned vector shape
    :param min_value: tuple of acceptable minimum value (size should be equal to shape param)
    :param max_value: tuple of acceptable maximum value (size should be equal to shape param)
    :param is_normalized: define if returned vector should be normalized
    :return: random vector
    """
    random_numbers = list()
    for pos in range(shape):
        random_numbers.append(random.uniform(min_value[pos], max_value[pos]))
    rand_vec = Vector2(*random_numbers)
    if is_normalized:
        vec_norm = norm(rand_vec)
        rand_vec = rand_vec / vec_norm
    return rand_vec


def decision(probability) -> bool:
    """
    Execute decision according probability

    :param probability: probability of successfully event
    :return: boolean indicating if event was a success
    """
    return random.uniform(0, 1) < probability


def random_number(min_value, max_value) -> float:
    """
    Calculate random number

    :param min_value: minimum acceptable value
    :param max_value: maximum acceptable value
    :return: random number
    """
    return random.uniform(min_value, max_value)
