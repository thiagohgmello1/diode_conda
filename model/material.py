import numpy as np
from scipy.constants import pi, h, epsilon_0, e


class Material:
    def __init__(
            self,
            mean_free_path: float,
            scalar_fermi_velocity: float,
            substrate_thickness: float,
            relax_time: float = None,
            gate_voltage: float = 10,
            permittivity: float = 1,
            permeability: float = 1,
            carrier_concentration=None
    ):
        self.relax_time = relax_time
        self.permittivity = permittivity
        self.permeability = permeability
        self.mean_free_path = mean_free_path
        self.scalar_fermi_velocity = scalar_fermi_velocity
        self.carrier_concentration = self._calc_carrier_concentration(
            gate_voltage, substrate_thickness, carrier_concentration
        )
        self.effective_mass = self._calc_effective_mass()
        self._calc_relax_time()


    def _calc_carrier_concentration(self, gate_voltage, substrate_thickness, carrier_concentration) -> float:
        """
        Define carrier concentration
        :param gate_voltage: voltage applied to the gate
        :param substrate_thickness: substrate thickness
        :param carrier_concentration: carrier concentration (if specified)
        :return: carrier concentration
        """
        if not carrier_concentration:
            carrier_concentration = epsilon_0 * self.permittivity * gate_voltage / (substrate_thickness * e)
        return carrier_concentration


    def _calc_effective_mass(self):
        """
        Calculate relative effective particle mass. From E=m_{e}c_{*}^2
        :return: effective particle mass
        """
        return (h * np.sqrt(self.carrier_concentration / pi)) / (2 * self.scalar_fermi_velocity)


    def _calc_relax_time(self):
        if not self.relax_time:
            self.relax_time = self.mean_free_path / self.scalar_fermi_velocity


if __name__ == '__main__':
    m = Material(1, 1e6, 300e-9, 10, 3.9)
    print(m.carrier_concentration)
