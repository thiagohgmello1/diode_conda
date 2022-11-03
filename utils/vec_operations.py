import numpy as np


def mirror(vec: np.array, normal: np.array):
    return vec - 2 * (np.dot(vec, normal)) / (np.linalg.norm(normal) ** 2) * normal


def calc_normal():
    pass
