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


def calc_asymmetry(current_list: list, volt: list) -> tuple[list, list]:
    """
    Calculate geometry asymmetry

    :param current_list: list of calculated currents
    :param volt: applied voltages
    :return: list of asymmetry by applied voltage
    """
    a = list()
    v = list()
    current_len = len(current_list)
    for pos in range(int(current_len / 2)):
        current_asymmetry = abs(current_list[pos] / current_list[current_len - pos - 1])
        v.append(volt[pos])
        a.append(current_asymmetry)

    return v, a


def plot_figs(asy_voltages, curr_voltages, asy, curr, drude_curr):
    fig_asymmetry = plt.figure(figsize=(12, 6))
    plt.plot(asy_voltages, asy)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.ylim(ymin=0)
    plt.savefig(f'outputs/asymmetry.png', dpi=fig_asymmetry.dpi)

    fig_iv = plt.figure(figsize=(12, 6))
    plt.plot(curr_voltages, curr, '--r', marker='o')
    if len(drude_curr) > 0:
        plt.plot(curr_voltages, drude_curr, '--', color='blue')
    plt.grid(True)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.savefig(f'outputs/currents.png', dpi=fig_iv.dpi)


def monte_carlo(voltage_list, topology, material, particle_model, geo, max_coll):
    for vol in voltage_list:
        voltages.append(-vol)
        vol = [vol, 0]
        e_field = Vector2(*vol) / (topology.bbox.xmax() - topology.bbox.xmin())
        system = System(
            topology=topology,
            material=material,
            particle=particle_model,
            electric_field=e_field,
            max_collisions=max_coll,
            max_time_simulation=material.relax_time
        )
        system.simulate(system.simulate_drude)
        simulation_current = system.cal_current()
        currents.append(simulation_current)
        save_current('outputs/currents.csv', simulation_current, geo, vol)

        print(f'Voltage: {-vol[0]}')
        print(f'Current:{simulation_current}')

        if 'rectangle' in geo:
            drude_current = drude_analytical_model(
                e_field=e_field,
                relax_time=material.relax_time,
                width=float(topology.bbox.ymax() - topology.bbox.ymin()),
                carrier_concentration=material.carrier_concentration,
                effective_mass=material.effective_mass * electron_mass
            )
            save_current('outputs/currents.csv', drude_current, geo, vol)
            drude_currents.append(drude_current)
            print(f'Drude current: {drude_current}')

        print(f'Time steps: {system.time_steps}')
        print(f'Collisions: {system.collisions}')
        print(f'-' * 100)
        print('\r')


if __name__ == '__main__':
    currents = list()
    drude_currents = list()
    voltages = list()

    MFPL = 200e-9
    mobility = 4
    carrier_c = None
    gate_voltage = 10
    f_velocity = c / 300
    sub_thickness = 300e-9
    material_permittivity = 3.9

    density = 50
    scale = 1e-7
    geometry = 'tests/diode12.svg'

    max_collisions = 100000
    voltage = np.linspace(-1, 1, num=11)

    mat = Material(
        mobility=mobility,
        mean_free_path=MFPL,
        gate_voltage=gate_voltage,
        carrier_concentration=carrier_c,
        scalar_fermi_velocity=f_velocity,
        substrate_thickness=sub_thickness,
        permittivity=material_permittivity
    )
    particle_m = Particle(density=density, effective_mass=mat.effective_mass, fermi_velocity=mat.scalar_fermi_velocity)
    pol = Topology.from_file(geometry, scale)

    exec_time = time.time()

    monte_carlo(voltage, pol, mat, particle_m, geometry, max_collisions)

    exec_time = time.time() - exec_time
    vol_asy, asymmetry = calc_asymmetry(currents, voltages)
    plot_figs(vol_asy, voltages, asymmetry, currents, drude_currents)
    print(f'Execution time: {exec_time}')
