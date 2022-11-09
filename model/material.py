import numpy as np
from scipy.constants import pi, h, electron_mass


class Material:
    def __init__(
            self,
            relax_time: float,
            scalar_fermi_velocity: float,
            carrier_concentration: float,
            permittivity: float = 1,
            permeability: float = 1
    ):
        self.relax_time = relax_time
        self.scalar_fermi_velocity = scalar_fermi_velocity
        self.effective_mass = self._calc_effective_mass(carrier_concentration)
        self.permittivity = permittivity
        self.permeability = permeability


    def _calc_effective_mass(self, carrier_concentration: float):
        return np.sqrt(h ** 2 * carrier_concentration / (4 * pi * electron_mass ** 2 * self.scalar_fermi_velocity ** 2))
