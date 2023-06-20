import json
import time

from pathlib import Path
from models.system import System
from solvers.mc_methods.drude import Drude
from solvers.monte_carlo import MonteCarlo


class Manager:
    def __init__(self, parse_args):
        self.parse_args = parse_args
        self.solver = None


    def set_method(self, system: System):
        """
        Set problem method. Ex.: drude
        :param system:
        :return: Chosen method
        """
        if self.parse_args.method == 'drude':
            return Drude(system)


    def set_solver(self, method, convergence_params):
        """
        Set problem solver. Ex.: mc (Monte Carlo)
        :param method: Previous chosen method
        :param convergence_params: convergence parameters
        :return:
        """
        if self.parse_args.solver == 'mc':
            self.solver = MonteCarlo(method, **convergence_params)


    def manage_files(self):
        """
        Manage selected files to be simulated
        :return:
        """
        total_time_exec = time.time()
        if self.parse_args.multi:
            self.multi_files()
        else:
            f = Path(f'parameters/{self.parse_args.single}.json')
            self.single_file(f)
        total_time_exec = time.time() - total_time_exec
        print(f'Total time exec: {"%s" % float("%.3g" % (total_time_exec / 60))} min')


    def single_file(self, file):
        """
        Simulate single file
        :param file:
        :return:
        """
        print(f'File: {file}')

        with open(file) as f:
            data = json.load(f)

        convergence_params = data['convergence']
        system = System(data['system'])
        method = self.set_method(system)
        self.set_solver(method, convergence_params)
        self.solver.create_voltage_range(**data['voltage'])
        voltages, currents, comp_currents, partial_exec_times = self.solver.simulate()
        self.save_data(voltages, currents, comp_currents, partial_exec_times, file.name)


    def multi_files(self):
        """
        Simulate multi files
        :return:
        """
        files = Path(f'parameters/{self.parse_args.multi}').glob('*')
        for file in files:
            print(f'File: {file}')
            self.single_file(file)


    def save_data(self, voltages: list, currents: list, comp_currents: list, partial_exec_times: list, file_name: str):
        """
        Save calculated current

        :param voltages: applied voltage
        :param currents: calculated current
        :param comp_currents: .svg simulated file
        :param partial_exec_times: id used to track simulation
        :param file_name:
        :return: None
        """
        strings_to_be_saved = list()
        for i in range(len(voltages)):
            if len(comp_currents) > 0:
                strings_to_be_saved.append(
                    f'{file_name};{currents[i]};{comp_currents[i]};{"%s" % float("%.3g" % voltages[i])};'
                    f'{self.parse_args.id};{partial_exec_times[i]}\n'
                )
            else:
                strings_to_be_saved.append(
                    f'{file_name};{currents[i]};{"%s" % float("%.3g" % voltages[i])};'
                    f'{self.parse_args.id};{partial_exec_times[i]}\n'
                )
        f = Path(f'outputs/{self.parse_args.output}.csv')
        if not f.is_file() and len(comp_currents) > 0:
            legend = 'file_name;sim_current;compare_current;voltage;id;exec_time\n'
            strings_to_be_saved.insert(0, legend)
        elif not f.is_file():
            legend = 'file_name;sim_current;voltage;id;exec_time\n'
            strings_to_be_saved.insert(0, legend)

        with open(f, 'a') as f:
            for row in strings_to_be_saved:
                f.write(row)

