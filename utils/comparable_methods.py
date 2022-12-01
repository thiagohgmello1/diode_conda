import numpy as np

from skgeom import Vector2
from scipy.constants import elementary_charge


def drude_analytical_model(
        width: float,
        relax_time: float,
        carrier_concentration: float,
        effective_mass: float,
        e_field: Vector2
):
    e_field = np.sqrt(float(e_field.squared_length()))
    return carrier_concentration * elementary_charge ** 2 * relax_time * e_field * width / effective_mass
