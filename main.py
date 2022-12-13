import time
import numpy as np
import matplotlib.pyplot as plt

from skgeom import Vector2
from scipy.constants import c, electron_mass
from model.particle import Particle
from model.topology import Topology
from model.material import Material
from simulators.monte_carlo import System, save_current
from utils.comparable_methods import drude_analytical_model


def calc_asymmetry(current_list: list, voltages: list) -> tuple[list, list]:
    """
    Calculate geometry asymmetry

    :param current_list: list of calculated currents
    :param voltages: applied voltages
    :return: list of asymmetry by applied voltage
    """
    asymmetry = list()
    vol = list()
    current_len = len(current_list)
    for pos in range(int(current_len / 2)):
        current_asymmetry = abs(current_list[current_len - pos - 1] / current_list[pos])
        vol.append(voltages[pos])
        asymmetry.append(current_asymmetry)

    return vol, asymmetry


if __name__ == '__main__':
    voltage = np.linspace(-0.5, 0.5, num=11)
    f_velocity = c / 300
    MFPL = 500e-9
    carrier_c = 7.2e15
    thickness = 300e-9
    gate_voltage = 10
    geometry = 'tests/diode7.svg'
    mat = Material(
        mean_free_path=MFPL,
        scalar_fermi_velocity=f_velocity,
        permittivity=3.9,
        substrate_thickness=thickness,
        gate_voltage=gate_voltage
    )
    particle_m = Particle(density=150, effective_mass=mat.effective_mass, fermi_velocity=mat.scalar_fermi_velocity)
    pol = Topology.from_file(geometry, 1e-7)
    currents = list()
    drude_currents = list()
    voltages = list()

    exec_time = time.time()

    for vol in voltage:
        voltages.append(-vol)
        vol = [vol, 0]
        e_field = Vector2(*vol) / (pol.bbox.xmax() - pol.bbox.xmin())
        system = System(
            particle=particle_m,
            topology=pol,
            material=mat,
            electric_field=e_field,
            max_collisions=100000,
            max_time_simulation=mat.relax_time
        )
        system.simulate(system.simulate_drude)
        simulation_current = system.cal_current()
        currents.append(simulation_current)
        save_current('outputs/currents.csv', simulation_current, geometry, vol)
        drude_current = drude_analytical_model(
            width=float(pol.bbox.ymax() - pol.bbox.ymin()),
            relax_time=mat.relax_time,
            carrier_concentration=mat.carrier_concentration,
            effective_mass=mat.effective_mass * electron_mass,
            e_field=e_field
        )
        save_current('outputs/currents.csv', drude_current, geometry, vol)
        drude_currents.append(drude_current)
        print(f'Current:{simulation_current}')
        print(f'Time steps: {system.time_steps}')
        print(f'Collisions: {system.collisions}')

    exec_time = time.time() - exec_time

    vol, asymmetry = calc_asymmetry(currents, voltages)
    fig_asymmetry = plt.figure(figsize=(12, 6))
    plt.plot(vol, asymmetry)
    plt.savefig('outputs/asymmetry.png', dpi=fig_asymmetry.dpi)

    fig_iv = plt.figure(figsize=(12, 6))
    plt.plot(voltages, currents)
    # plt.plot(voltages, drude_currents, color='blue')
    plt.savefig('outputs/currents.png', dpi=fig_iv.dpi)

    print(f'Execution time: {exec_time}')
