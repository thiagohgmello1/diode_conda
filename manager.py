import json
import time

from pathlib import Path
from datetime import datetime
from models.system import System
from solvers.mc_methods.drude import Drude
from solvers.monte_carlo import MonteCarlo


class Manager:
    def __init__(self, parse_args):
        self.parse_args = parse_args
        self.solver = None


    def set_method(self, system: System):
        if self.parse_args.method == 'drude':
            return Drude(system)


    def set_solver(self, method, convergence_params):
        if self.parse_args.solver == 'mc':
            self.solver = MonteCarlo(method, **convergence_params)


    def manage_files(self):
        if self.parse_args.multi:
            self.multi_files()
        else:
            f = Path(f'parameters/{self.parse_args.single}.json')
            self.single_file(f)


    def single_file(self, file):
        print(f'File: {file}')

        with open(file) as f:
            data = json.load(f)

        convergence_params = data['convergence']
        total_time_exec = time.time()
        system = System(data['system'])
        method = self.set_method(system)
        self.set_solver(method, convergence_params)
        self.solver.create_voltage_range(**data['voltage'])
        voltages, currents, comp_currents = self.solver.simulate()
        print(f'Total time exec: {"%s" % float("%.3g" % (total_time_exec / 60))} min')
        # self.save_current(voltages, currents, comp_currents, id_tracker)


    def multi_files(self):
        files = Path(f'parameters/{self.parse_args.multi}').glob('*')
        total_time_exec = time.time()
        for file in files:
            print(f'File: {file}')
            self.single_file(file)
        total_time_exec = time.time() - total_time_exec
        print(f'Total time exec: {"%s" % float("%.3g" % (total_time_exec / 60))} min')


    def save_current(self, meas_current: float, simulated_geometry: str, volt: float, id_tracker: str):
        """
        Save calculated current

        :param volt: applied voltage
        :param meas_current: calculated current
        :param simulated_geometry: .svg simulated file
        :param id_tracker: id used to track simulation
        :return: None
        """
        now = datetime.now()
        date_string = now.strftime("%d/%m/%Y %H:%M:%S")
        with open(self.parse_args.output, 'a') as f:
            string_to_be_saved = \
                f'{date_string};{simulated_geometry};{meas_current};{"%s" % float("%.1g" % volt)};{id_tracker}\n'
            f.write(string_to_be_saved)

