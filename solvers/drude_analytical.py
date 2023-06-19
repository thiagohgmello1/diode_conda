import numpy as np

from scipy.constants import elementary_charge
from models.electric_field import ElectricField


def drude_analytical_model(
        width: float,
        relax_time: float,
        carrier_concentration: float,
        effective_mass: float,
        e_field: ElectricField
):
    e_field = -1 * np.sign(float(e_field.vector([0, 0]).x())) * np.sqrt(float(e_field.vector([0, 0]).squared_length()))
    return carrier_concentration * elementary_charge ** 2 * relax_time * e_field * width / effective_mass
