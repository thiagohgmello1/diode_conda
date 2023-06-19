import argparse

from manager import Manager


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for diode Monte Carlo simulation')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--multi', '--m', type=str, help="Define simulation directory for multi-files")
    group.add_argument('--single', '--s', type=str, help="Define simulation file for single-simulation")
    parser.add_argument('--output', '--o', default='currents', type=str, help='Current save output file name')
    parser.add_argument('--id', type=str, help='ID code used to track simulation. Can be a string without whitespace')
    parser.add_argument('--method', default='drude', type=str, help='Numerical method')
    parser.add_argument('--solver', default='mc', type=str, help='solver')

    args = parser.parse_args()

    exec_time_list = list()
    asymmetry_list = list()
    curr_list = list()
    drude_curr_list = list()
    vol_asy_list = list()
    voltages_list = list()

    manager = Manager(args)
    manager.manage_files()
