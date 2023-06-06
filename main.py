import time
import json
import argparse
import numpy as np

from pathlib import Path
from model.particle import Particle
from model.topology import Topology
from model.material import Material
from simulators.monte_carlo import monte_carlo
from utils.post_processing import calc_asymmetry, plot_figs


def create_voltage_range(v_min, v_max, num_points):
    if num_points == 1:
        return [v_max]
    else:
        return np.linspace(v_min, v_max, num=num_points)


def chose_topology(geometry_dict) -> Topology:
    geometry_type = geometry_dict['input_style']
    del geometry_dict['input_style']
    if geometry_type == 'file':
        return Topology.from_file(**geometry_dict)
    else:
        return Topology.from_points(**geometry_dict)


def simulate(file_name: Path, out_file, id_tracker):
    currents = list()
    drude_currents = list()
    voltages = list()

    with open(file_name) as f:
        data = json.load(f)

    voltage = create_voltage_range(**data['voltage'])

    mat = Material(**data['material'])
    particle_m = Particle(
        **data['particle'], effective_mass=mat.effective_mass, fermi_velocity=mat.scalar_fermi_velocity
    )
    pol = chose_topology(data['geometry'])

    exec_time = time.time()

    monte_carlo(voltage, pol, mat, particle_m, voltages, currents, drude_currents, **data['convergence'],
                out_file=out_file, id_tracker=id_tracker)
    exec_time = time.time() - exec_time

    vol_asy, asymmetry = calc_asymmetry(currents, voltages)

    return exec_time, vol_asy, asymmetry, voltages, currents, drude_currents


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for diode Monte Carlo simulation')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--multi', '--m', type=str, help="Define simulation directory for multi-files")
    group.add_argument('--single', '--s', type=str, help="Define simulation file for single-simulation")
    parser.add_argument('--output', '--o', default='currents', type=str, help='Current save output file name')
    parser.add_argument('--id', type=str, help='ID code used to track simulation. Can be a string without whitespace')

    args = parser.parse_args()

    exec_time_list = list()
    asymmetry_list = list()
    curr_list = list()
    drude_curr_list = list()
    vol_asy_list = list()
    voltages_list = list()

    if args.multi:
        files = Path(f'parameters/{args.multi}').glob('*')
        for file_sim in files:
            exec_time_aux, vol_asy_aux, asymmetry_aux, voltages_aux, curr_aux, drude_curr_aux = \
                simulate(file_sim, args.output, args.id)
            exec_time_list.append(exec_time_aux)
            vol_asy_list.append(vol_asy_aux)
            asymmetry_list.append(asymmetry_aux)
            voltages_list.append(voltages_aux)
            curr_list.append(curr_aux)
            drude_curr_list.append(drude_curr_aux)
    else:
        file = Path(f'parameters/{args.single}.json')
        exec_time_aux, vol_asy_aux, asymmetry_aux, voltages_aux, curr_aux, drude_curr_aux = \
            simulate(file, args.output, args.id)
        plot_figs(vol_asy_aux, voltages_aux, asymmetry_aux, curr_aux, drude_curr_aux)
        print(f'Execution time: {"%s" % float("%.3g" % (exec_time_aux / 60))} min')
