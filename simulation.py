#! /usr/bin/env python

# ===============================
# >>> IMPORTS
# ===============================

import time
import math

import pysic
import friction_tools as ft


# ===============================
# >>> CONFIGS
# ===============================

# ===============================
# >>> FUNCTIONS
# ===============================


def create_lj_interactions(strengths):
    """
    Creates Lernnard-Jones potential interactions between materials.
    @param strengths: {'material': 1.0} - dict of material-strength pairs.
    """
    sim.create_interaction(
        [fixed_slab_material, fixed_slab_material],
        strengths.get(fixed_slab_material, 0.0),
        10 / (fixed_slab_xy * (2 ** 0.5))
    )
    sim.create_interaction(
        [moving_slab_material, moving_slab_material],
        strengths.get(moving_slab_material, 0.0),
        10 / (moving_slab_xy * (2 ** 0.5))
    )
    sim.create_interaction(
        [fixed_slab_material, moving_slab_material],
        1.0,
        2.0
    )


def create_coulomb_interactions(charges):
    """
    Creates Coulombs summation interactions between materials.
    @param charges: {'material': 1.0} - dict of material-charge pairs.
    """
    def charge_by_elements(element):
        """
        Maps elements with their according charges.
        """
        return charges.get(element, 0.0)

    # Set charges
    sim.system.set_charges(
        map(charge_by_elements, sim.system.get_chemical_symbols())
    )
    # Create instance of the CoulombSummation class
    coulomb_sum = pysic.CoulombSummation()
    # Find good enough computational parameters
    parameters = pysic.interactions.coulomb.estimate_ewald_parameters()
    # Set the parameters to the summation object
    coulomb_sum.set_parameter_values(list(parameters))
    # Add the created CoulombSummation object to the Pysic calculator object
    sim.calc.set_coulomb_summation(coulomb_sum)


def friction_initial_velocity_x(threshold_x, initial_velocity_x=0.1):
    """
    Measures friction by giving the moving slab an initial velocity in
    x-direction and simulating until the slab has stopped moving.
    @param threshold_x: 1.0 - values below threshold count as being stopped.
    """
    # Set initial velocity in x-direction
    sim.set_velocities(
        indices=moving_mat_indices,
        velocity=[initial_velocity_x, 0, 0]
    )
    print('Slab is now moving. Simming until it stops...')
    while sim.get_avg_velocity(moving_mat_indices)[0] > threshold_x:
        print('Velocities are now: '+str(
            sim.get_avg_velocity(moving_mat_indices)
        ))
        sim.run_simulation(300)

    print('Velocities are now: '+str(
        sim.get_avg_velocity(moving_mat_indices)
    ))


def friction_constant_velocity_x(indices, constant_velocity_x):

    sim.fix_velocities(
        indices=indices,
        velocity=[constant_velocity_x, 0, 0]
    )

    print('Slab is now moving in constant velocity...')

    sim.run_simulation(1000)


# ===============================
# >>> MAIN
# ===============================

def main(strengths):
    moving_slab_material = 'Ag'
    moving_slab_xy = 3  # 3
    moving_slab_z = 2  # 2

    fixed_slab_material = 'Au'
    fixed_slab_xy = 5  # 5
    fixed_slab_z = 2  # 2

    element_charges = {
        'Au': 0.01,
        'Ag': 0.005
    }

    element_lj_strengths = strengths
    """= {
        'Au': 1.0,
        'Ag': 0.5
    }"""

    # Create simulation object
    sim = ft.FrictionSimulation()


    # Calculate thresholds for being in touch with the surface
    contact_threshold = math.ceil(10.0 / float(moving_slab_xy)) * moving_slab_z
    velocity_threshold_x = 0.04


    # Create the two slabs and store their indices
    # TODO: Edit to use z-config
    sim.create_slab(
        element=fixed_slab_material,
        xy_cells=fixed_slab_xy,
        z_cells=2,
        top_z=0.0
    )
    sim.create_slab(
        element=moving_slab_material,
        xy_cells=moving_slab_xy,
        z_cells=2,
        bottom_z=5.0
    )

    fixed_mat_indices = sim.get_indices_by_element(fixed_slab_material)
    moving_mat_indices = sim.get_indices_by_element(moving_slab_material)

    # TODO: Edit to use z-config
    bottom_indices = sim.get_indices_z_less_than(-2)  # Bottom row of fixed materials atoms
    top_indices = sim.get_indices_z_more_than(8.0)  # Top row of moving materials atoms


    # Create interactions
    create_lj_interactions(element_lj_strengths)
    # create_coulomb_interactions(element_charges)  # Currently too slow!


    # Create dynamics and configure the setup
    sim.create_dynamics(dt=1.0, temperature=300)

    # Fix the bottom slab and make the top slab move down
    sim.fix_positions(bottom_indices)
    sim.add_constant_force(top_indices, [0, 0, -0.1])
    # sim.fix_velocities(top_indices, [0, 0, -0.005], [True, True, True])

    sim.set_temperature(300)  # Try to reset the temperature


    # Save tmps, print stats and time the simulation
    sim.save_trajectory_during_simulation(interval=50.0)
    sim.print_stats_during_simulation(interval=50.0)
    t0 = time.time()

    print('Starting the simulation')
    print(sim.get_avg_forces(moving_mat_indices))
    while len(sim.get_indices_z_more_than(contact_threshold)):
        print(len(sim.get_indices_z_more_than(contact_threshold)))
        print(sim.get_avg_velocity(moving_mat_indices))
        print(sim.get_avg_forces(moving_mat_indices))
        sim.run_simulation(500)

    print('Moving slab hit the fixed slab.')
    sim.remove_constraints()
    print('Starting friction measurement phase.')
    xpos0 = [p[0] for p in sim.get_positions(indices=moving_mat_indices)]
    friction_initial_velocity_x(0.005, 0.1)
    xpos1 = [p[0] for p in sim.get_positions(indices=moving_mat_indices)]
    xdiff = []
    for i in range(0, len(xpos0)):
        xdiff.append(xpos1[i] - xpos0[i])

    out = open('ownpos.txt', 'w')
    for i in range(0, len(xpos0)):
        out.write('{0}\n'.format(xpos0[i], xpos1[i]))
    out.close()

    print('All done, generating simulation file.')
    print(sim.get_avg_forces(moving_mat_indices))
    sim.write_positions_to_file(indices=moving_mat_indices)
    sim.write_average_position_to_file(indices=moving_mat_indices)
    sim.write_velocities_to_file(indices=moving_mat_indices)
    ft.trajectory_to_xyz()


#=====================
#>>>>> MAIN 2 <<<<<<<<
#=====================


#sim = ft.FrictionSimulation()


## Calculate thresholds for being in touch with the surface
#contact_threshold = math.ceil(10.0 / float(moving_slab_xy)) * moving_slab_z
#velocity_threshold_x = 0.04


## Create the two slabs and store their indices
## TODO: Edit to use z-config
#sim.create_slab(element=fixed_slab_material, xy_cells=fixed_slab_xy, z_cells=2, top_z=0.0)
#sim.create_slab(element=moving_slab_material, xy_cells=moving_slab_xy, z_cells=2, bottom_z=5.0)

#fixed_mat_indices = sim.get_indices_by_element(fixed_slab_material)
#moving_mat_indices = sim.get_indices_by_element(moving_slab_material)

## TODO: Edit to use z-config
#bottom_indices = sim.get_indices_z_less_than(-2)  # Bottom row of fixed materials atoms
#top_indices = sim.get_indices_z_more_than(8.0)  # Top row of moving materials atoms


## Create interactions
#create_lj_interactions(element_lj_strengths)
#create_coulomb_interactions(element_charges)  # Currently too slow!


## Create dynamics and configure the setup
#sim.create_dynamics(dt=1.0, temperature=300)

## Fix the bottom slab and make the top slab move down
#sim.fix_positions(bottom_indices)
#sim.add_constant_force(top_indices, [0, 0, -0.1])
## sim.fix_velocities(top_indices, [0, 0, -0.005], [True, True, True])

#sim.set_temperature(300)  # Try to reset the temperature


## Save tmps, print stats and time the simulation
#sim.save_trajectory_during_simulation(interval=50.0)
#sim.print_stats_during_simulation(interval=50.0)
#t0 = time.time()

#print('Starting the simulation')

#while len(sim.get_indices_z_more_than(contact_threshold)):
    #print(len(sim.get_indices_z_more_than(contact_threshold)))
    #print(sim.get_avg_velocity(moving_mat_indices))
    #sim.run_simulation(500)

#print('Moving slab hit the fixed slab.')

#print('Starting friction measurement phase.')
#friction_constant_velocity_x([initial_velocity, 0.0, 0.0])

#print('All done, generating simulation file.')
#ft.trajectory_to_xyz()