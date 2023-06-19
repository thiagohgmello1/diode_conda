import matplotlib.pyplot as plt

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
        current_asymmetry = abs(current_list[pos] / current_list[current_len - pos - 1])
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


def plot_stable_current(currents, voltage):
    fig_curr, ax = plt.subplots(figsize=(12, 6))
    formatter = EngFormatter(places=1, sep="\N{THIN SPACE}")
    ax.yaxis.set_major_formatter(formatter)
    plt.plot(currents)
    ax.set_xlabel('time step')
    ax.set_ylabel('Current [A]')
    plt.savefig(f"outputs/current_stable/currents{'%s' % float('%.1g' % voltage)}.png", dpi=fig_curr.dpi)
