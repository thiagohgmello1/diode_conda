import time
import json
import numpy as np

from model.particle import Particle
from model.topology import Topology
from model.material import Material
from simulators.monte_carlo import monte_carlo
from utils.post_processing import calc_asymmetry, plot_figs


def create_voltage_range(v_min, v_max, num_points):
    return np.linspace(v_min, v_max, num=num_points)


def chose_topology(geometry_dict) -> Topology:
    geometry_type = geometry_dict['input_style']
    del geometry_dict['input_style']
    if geometry_type == 'file':
        return Topology.from_file(**geometry_dict)
    else:
        return Topology.from_points(**geometry_dict)


if __name__ == '__main__':
    currents = list()
    drude_currents = list()
    voltages = list()

    with open('parameters/rectangle.json') as f:
        data = json.load(f)

    voltage = create_voltage_range(**data['voltage'])

    mat = Material(**data['material'])
    particle_m = Particle(
        **data['particle'], effective_mass=mat.effective_mass, fermi_velocity=mat.scalar_fermi_velocity
    )
    pol = chose_topology(data['geometry'])

    exec_time = time.time()

    monte_carlo(voltage, pol, mat, particle_m, voltages, currents, drude_currents, **data['convergence'])
    exec_time = time.time() - exec_time

    vol_asy, asymmetry = calc_asymmetry(currents, voltages)
    plot_figs(vol_asy, voltages, asymmetry, currents, drude_currents)
    print(f'Execution time: {exec_time}')
