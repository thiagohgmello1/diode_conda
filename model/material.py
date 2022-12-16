import numpy as np
from scipy.constants import pi, h, epsilon_0, elementary_charge, electron_mass


class Material:
    def __init__(
            self,
            mean_free_path: float,
            scalar_fermi_velocity: float,
            mobility: float = None,
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
        :param mobility: electron mobility [m^2/(Vs)]
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
        self.mobility = self._calc_mobility(mobility)

    def _calc_carrier_concentration(self, gate_voltage, substrate_thickness, carrier_concentration) -> float:
        """
        Define carrier concentration
        :param gate_voltage: voltage applied to the gate
        :param substrate_thickness: substrate thickness
        :param carrier_concentration: carrier concentration (if specified)
        :return: carrier concentration
        """
        if not carrier_concentration:
            carrier_concentration = epsilon_0 * self.permittivity * gate_voltage / \
                                    (substrate_thickness * elementary_charge)
        return carrier_concentration

    def _calc_effective_mass(self):
        """
        Calculate relative effective particle mass. From E=m_{e}v_{*}^2

        :return: effective particle mass
        """
        return (h * np.sqrt(self.carrier_concentration / pi)) / (2 * self.scalar_fermi_velocity * electron_mass)

    def _calc_relax_time(self):
        """
        Calculate relaxation time (if not specified)

        :return: relaxation time
        """
        if not self.relax_time:
            self.relax_time = self.mean_free_path / self.scalar_fermi_velocity


    def _calc_mobility(self, mobility):
        """
        Calculate material mobility

        :param mobility: input mobility (if specified)
        :return: mobility
        """
        if not mobility:
            return elementary_charge * self.relax_time / (self.effective_mass * electron_mass)
        else:
            return mobility
