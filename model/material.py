import numpy as np
from scipy.constants import pi, h, hbar, k, epsilon_0, elementary_charge, electron_mass


class Material:
    def __init__(
            self,
            carrier_concentration: dict,
            scalar_fermi_velocity: float,
            mean_free_path: float,
            mobility: float = None,
            relax_time: float = None,
            permittivity: float = 1,
            permeability: float = 1,
            effective_mass: float = None
    ):
        """
        Represent the chosen material. All values must be in international system (m, s, etc.)

        :param carrier_concentration: dict for select carrier concentration calculu method
        :param scalar_fermi_velocity: defined scalar Fermi velocity (m/s)
        :param mobility: electron mobility [m^2/(Vs)]
        :param mean_free_path: material mean free path (Lambda_{MFP})
        :param relax_time: relaxation time (s) (if desired to define)
        :param permittivity: electric permittivity (F/m)
        :param permeability: magnetic permeability (H/m)
        """
        self.scalar_fermi_velocity = scalar_fermi_velocity
        self.mean_free_path = mean_free_path
        self.relax_time = self._calc_relax_time(relax_time)
        self.permittivity = permittivity
        self.permeability = permeability
        self.carrier_concentration = self._calc_carrier_concentration(carrier_concentration)
        self.effective_mass = self._calc_effective_mass(effective_mass)
        self.mobility = self._calc_mobility(mobility)


    def _calc_carrier_concentration(self, carrier_dict) -> float:
        """
        Define carrier concentration
        :param carrier_dict: dictionary with carrier concentration data
        :return: carrier concentration
        """
        if carrier_dict['method'] == 'gate_voltage':
            carrier_concentration = epsilon_0 * self.permittivity * carrier_dict['gate_voltage'] / \
                                    (carrier_dict['substrate_thickness'] * elementary_charge)
        elif carrier_dict['method'] == 'temp':
            carrier_concentration = pi / 6 * (k * carrier_dict['temp'] / (hbar * self.scalar_fermi_velocity)) ** 2
        else:
            carrier_concentration = carrier_dict['value']

        return carrier_concentration


    def _calc_effective_mass(self, effective_mass):
        """
        Calculate relative effective particle mass. From E=m_{e}v_{*}^2

        :return: effective particle mass
        """
        if not effective_mass:
            return (h * np.sqrt(self.carrier_concentration / pi)) /\
                                  (2 * self.scalar_fermi_velocity * electron_mass)
        else:
            return effective_mass


    def _calc_mean_free_path(self, mean_free_path):
        if not mean_free_path:
            self.mean_free_path = self.relax_time * self.scalar_fermi_velocity
        else:
            self.mean_free_path = mean_free_path


    def _calc_relax_time(self, relax_time):
        """
        Calculate relaxation time (if not specified)

        :return: relaxation time
        """
        if not relax_time:
            return self.mean_free_path / self.scalar_fermi_velocity
        else:
            return relax_time


    def _calc_mobility(self, mobility):
        """
        Calculate material mobility

        :param mobility: input mobility (if specified)
        :return: None
        """
        if not mobility:
            return elementary_charge * self.relax_time / (self.effective_mass * electron_mass)
        else:
            return mobility
