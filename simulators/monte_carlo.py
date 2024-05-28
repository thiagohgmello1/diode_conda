from skgeom import Vector2
from model.system import System
from scipy.constants import electron_mass
from utils.post_processing import save_current
from matplotlib.ticker import EngFormatter
from simulators.drude_analytical import drude_analytical_model


def monte_carlo(
        volt,
        topology,
        material,
        particle_model,
        max_coll,
        n_particles=100,
        check_condition='time',
        plot_current=True
):
    volt_vec = [-volt, 0]
    # For now, simulator considers only x electric fields
    e_field = Vector2(*volt_vec) / (topology.bbox.xmax() - topology.bbox.xmin())
    system = System(
        topology=topology,
        material=material,
        particle=particle_model,
        electric_field=e_field,
        max_collisions=max_coll,
        number_of_particles=n_particles,
        check_condition=check_condition
    )
    system.simulate(system.simulate_drude, volt, plot_current)
    simulation_current = system.cal_current()
    return e_field, simulation_current, system.time_steps_count, system.collisions_count


def monte_carlo_non_opt(
        voltage_range,
        topology,
        material,
        particle_model,
        voltages,
        currents,
        drude_currents,
        geo,
        max_coll,
        n_particles=100,
        out_file='currents',
        id_tracker='test',
        check_condition='time'
):
    for volt in voltage_range:
        eng_formatter = EngFormatter(places=4, unit='A')
        voltages.append(volt)
        e_field, simulation_current, time_steps_count, collisions_count = monte_carlo(
            volt, topology, material, particle_model, max_coll, n_particles, check_condition
        )
        currents.append(simulation_current)
        save_current(f'outputs/{out_file}.csv', simulation_current, geo, volt, id_tracker)

        print(f"Voltage: {'%s' % float('%.1g' % volt)}")
        print(f"Current:{eng_formatter.format_eng(num=simulation_current)}A")

        if 'rectangle' in geo:
            drude_current = drude_analytical_model(
                e_field=e_field,
                relax_time=material.relax_time,
                width=float(topology.bbox.ymax() - topology.bbox.ymin()),
                carrier_concentration=material.carrier_concentration,
                effective_mass=material.effective_mass * electron_mass
            )
            save_current(f'outputs/{out_file}.csv', drude_current, geo, volt, id_tracker)
            drude_currents.append(drude_current)
            print(f'Drude current: {drude_current}')

        print(f'Time steps: {time_steps_count}')
        print(f'Collisions: {collisions_count}')
        print(f'-' * 100)
        print('\r')
