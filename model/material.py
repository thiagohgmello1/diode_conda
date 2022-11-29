import numpy as np
from scipy.constants import pi, h, epsilon_0, e, c


class Material:
    def __init__(
            self,
            mean_free_path: float,
            scalar_fermi_velocity: float,
            substrate_thickness: float = 1,
            relax_time: float = None,
            gate_voltage: float = 10,
            permittivity: float = 1,
            permeability: float = 1,
            carrier_concentration=None
    ):
        """
        Class responsible to represent used material. All values must be in international system (m, s, etc.)

        :param mean_free_path: material mean free path (Lambda_{MFP})
        :param scalar_fermi_velocity: defined scalar Fermi velocity
        :param substrate_thickness: substrate thickness (responsible from carrier concentration)
        :param relax_time: relaxation time (if desired to define)
        :param gate_voltage: applied gate voltage
        :param permittivity: electric permittivity
        :param permeability: magnetic permeability
        :param carrier_concentration: material carrier concentration (if desired to define)
        """
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
        Calculate relative effective particle mass. From E=m_{e}v_{*}^2

        :return: effective particle mass
        """
        return (h * np.sqrt(self.carrier_concentration / pi)) / (2 * self.scalar_fermi_velocity)


    def _calc_relax_time(self):
        """
        Calculate relaxation time (if not specified)

        :return: relaxation time
        """
        if not self.relax_time:
            self.relax_time = self.mean_free_path / self.scalar_fermi_velocity


if __name__ == '__main__':
    f_velocity = c / 300
    mean_FPL = 200e-9
    carrier_c = 1.1e16
    m = Material(mean_free_path=mean_FPL, scalar_fermi_velocity=f_velocity, carrier_concentration=carrier_c)
    print(m.carrier_concentration)
