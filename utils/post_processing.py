import matplotlib.pyplot as plt

from datetime import datetime
from matplotlib.ticker import EngFormatter


def calc_asymmetry(current_list: list, voltages: list) -> tuple[list, list]:
    """
    Calculate geometry asymmetry

    :param current_list: list of calculated currents
    :param voltages: applied voltages
    :return: list of asymmetry by applied voltage
    """
    a = list()
    v = list()
    current_len = len(current_list)
    for pos in range(int(current_len / 2)):
        current_asymmetry = abs(current_list[current_len - pos - 1] / current_list[pos])
        v.append(voltages[pos])
        a.append(current_asymmetry)

    return v, a


def plot_figs(asymmetric_voltages, curr_voltages, asymmetry, current, drude_curr):
    fig_asymmetry = plt.figure(figsize=(12, 6))
    plt.plot(asymmetric_voltages, asymmetry)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.ylim(ymin=0)
    plt.savefig(f'outputs/asymmetry.png', dpi=fig_asymmetry.dpi)

    fig_iv = plt.figure(figsize=(12, 6))
    plt.plot(curr_voltages, current, '--r', marker='o')
    if len(drude_curr) > 0:
        plt.plot(curr_voltages, drude_curr, '--', color='blue')
    plt.grid(True)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.savefig(f'outputs/currents.png', dpi=fig_iv.dpi)


def save_current(current_file: str, meas_current: float, simulated_geometry: str, volt: float, id_tracker: str):
    """
    Save calculated current

    :param volt: applied voltage
    :param current_file: output file
    :param meas_current: calculated current
    :param simulated_geometry: .svg simulated file
    :param id_tracker: id used to track simulation
    :return: None
    """
    now = datetime.now()
    date_string = now.strftime("%d/%m/%Y %H:%M:%S")
    with open(current_file, 'a') as f:
        string_to_be_saved = \
            f'{date_string};{simulated_geometry};{meas_current};{"%s" % float("%.1g" % volt)};{id_tracker}\n'
        f.write(string_to_be_saved)


def progress_bar(progress, total):
    """
    Print progress bar

    :param progress: evolved situation
    :param total: expect max situation
    :return: None
    """
    percent = 100 * (progress / total)
    bar = '█' * int(percent) + '-' * (100 - int(percent))
    print(f'\r|{bar}| {percent:.2f}%', end='\r')


def plot_stable_current(currents, voltage):
    fig_curr, ax = plt.subplots(figsize=(12, 6))
    formatter = EngFormatter(places=1, sep="\N{THIN SPACE}")
    ax.yaxis.set_major_formatter(formatter)
    plt.plot(currents)
    ax.set_xlabel('time step')
    ax.set_ylabel('Current [A]')
    plt.savefig(f"outputs/current_stable/currents{'%s' % float('%.1g' % voltage)}.png", dpi=fig_curr.dpi)
    plt.close(fig_curr)
