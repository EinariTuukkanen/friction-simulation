#! /usr/bin/env python
# coding=utf-8
"""A module containing routines for setting up simple friction simulations.
"""

# get Pysic - this is the package that handles atomic interactions
# see http://thynnine.github.io/pysic/
import pysic

# get some components of ASE - this is the package that handles dynamics
# see https://wiki.fysik.dtu.dk/ase/
from ase.lattice.compounds import Rocksalt
#from ase.structure import bulk
from ase.lattice import bulk
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from ase.md.nvtberendsen import NVTBerendsen
from ase.optimize import BFGSLineSearch
from ase.io.trajectory import PickleTrajectory
from ase.constraints import FixAtoms
import ase.io
from ase import Atoms
from ase import units

# get numeric tools from NumPy - see http://scipy.org
import numpy as np
import numpy.random as nr

# get Python tools - these are part of the Python default library
import math
import time
import copy


# get plotting tools
try:
    import matplotlib.pyplot as plt
except:
    print "MatPlotLib not found, plotting will not work"


def set_random_number_seed(seed=5801):
    """Initialize random number generator (RNG).

    By default, the RNG is always given the same seed so that simulations are reproducible.
    However, if you would like to rerun simulations with different random sequences,
    give a new random seed here.

    This affects :meth:`set_initial_temperature` and :meth:`create_dynamics` / :math:`run_simulation`
    if Langevin dynamics are chosen.

    Parameters:

    seed: integer
        a new RNG seed
    """
    nr.seed(seed)

set_random_number_seed()

# mod the Langevin dynamics class from ASE
from ase.parallel import world 
from ase.md.md import MolecularDynamics 
import sys

class LangevinMod(Langevin): 
    
    def __init__(self, atoms, timestep, temperature, friction, fixcm=False, 
                   trajectory=None, logfile=None, loginterval=1, communicator=world,
                   applyArray=None):
          """A copy of the Langevin dynamics class from ASE in case modifications need to be made.

          In Langevin dynamics, a friction force is applied to moving atoms to prevent
          them from accumulating too much kinetic energy. This will drive the temperature
          to zero over time, and so it is balanced by a randomly fluctuating force that
          affects each atom independently. The strength of this random force will determine
          the temperature of the system in equilibrium.
          """
          MolecularDynamics.__init__(self, atoms, timestep, trajectory, 
                                logfile, loginterval) 
          self.temp = temperature
          self.frict = friction 
          self.fixcm = fixcm  # will the center of mass be held fixed? 
          self.communicator = communicator 
          self.updatevars() 
          self.apply = applyArray
          
    def updatevars(self): 
          dt = self.dt
          # If the friction is an array some other constants must be arrays too. 
          self._localfrict = hasattr(self.frict, 'shape') 
          lt = self.frict * dt 
          masses = self.masses 
          sdpos = dt * np.sqrt(self.temp / masses.reshape(-1) * (2.0/3.0 - 0.5 * lt) * lt) 
          sdpos.shape = (-1, 1) 
          sdmom = np.sqrt(self.temp * masses.reshape(-1) * 2.0 * (1.0 - lt) * lt) 
          sdmom.shape = (-1, 1) 
          pmcor = np.sqrt(3.0)/2.0 * (1.0 - 0.125 * lt) 
          cnst = np.sqrt((1.0 - pmcor) * (1.0 + pmcor)) 
   
          act0 = 1.0 - lt + 0.5 * lt * lt 
          act1 = (1.0 - 0.5 * lt + (1.0/6.0) * lt * lt) 
          act2 = 0.5 - (1.0/6.0) * lt + (1.0/24.0) * lt * lt 
          c1 = act1 * dt / masses.reshape(-1) 
          c1.shape = (-1, 1) 
          c2 = act2 * dt * dt / masses.reshape(-1) 
          c2.shape = (-1, 1) 
          c3 = (act1 - act2) * dt 
          c4 = act2 * dt 
          del act1, act2 
          if self._localfrict: 
              # If the friction is an array, so are these 
              act0.shape = (-1, 1) 
              c3.shape = (-1, 1) 
              c4.shape = (-1, 1) 
              pmcor.shape = (-1, 1) 
              cnst.shape = (-1, 1) 
          self.sdpos = sdpos 
          self.sdmom = sdmom 
          self.c1 = c1 
          self.c2 = c2 
          self.act0 = act0 
          self.c3 = c3 
          self.c4 = c4 
          self.pmcor = pmcor 
          self.cnst = cnst 
          self.natoms = self.atoms.get_number_of_atoms() 
          
          
    def step(self, f): 
          atoms = self.atoms 
          p = self.atoms.get_momenta() 
   
          random1 = nr.standard_normal(size=(len(atoms), 3)) 
          random2 = nr.standard_normal(size=(len(atoms), 3)) 
   
          if self.communicator is not None: 
              self.communicator.broadcast(random1, 0) 
              self.communicator.broadcast(random2, 0) 
           
          rrnd = self.sdpos * random1 
          prnd = (self.sdmom * self.pmcor * random1 + 
                  self.sdmom * self.cnst * random2) 
   
          if self.fixcm: 
              rrnd = rrnd - np.sum(rrnd, 0) / len(atoms) 
              prnd = prnd - np.sum(prnd, 0) / len(atoms) 
              rrnd *= np.sqrt(self.natoms / (self.natoms - 1.0)) 
              prnd *= np.sqrt(self.natoms / (self.natoms - 1.0)) 
   
          atoms.set_positions(atoms.get_positions() + 
                              self.c1 * p + 
                              self.c2 * f + rrnd) 
          p *= self.act0 
          p += self.c3 * f + prnd 
          atoms.set_momenta(p) 
                         
          f = atoms.get_forces() 
          atoms.set_momenta(atoms.get_momenta() + self.c4 * f) 
          return f 
  
  
  
        

class FixPositions:
    """A constraint class for fixing atomic positions
    """
    def __init__(self, indices, fixed):
        self.indices = indices
        self.fixed = fixed

    def adjust_positions(self, oldpositions, newpositions):
        for index in self.indices:
            for i in range(3):
                if self.fixed[i]:
                    newpositions[index,i] = oldpositions[index,i]

    def adjust_forces(self, positions, forces):
        for index in self.indices:
            for i in range(3):
                if self.fixed[i]:
                    forces[index,i] = 0.0

    def copy(self):
        return copy.deepcopy(self)


class FixVelocities:
    """A constraint class for fixing atomic velocities
    """
    def __init__(self, indices, velocity, dt, fixed):
        self.indices = indices
        self.fixed = fixed
        self.step = np.array(velocity, dtype=float)*dt
        self.velocity = np.array(velocity)

    def adjust_positions(self, oldpositions, newpositions):
        for index in self.indices:
            for i in range(3):
                if self.fixed[i]:
                    newpositions[index,i] = oldpositions[index,i] + self.step[i]

    def adjust_forces(self, positions, forces):
        for index in self.indices:
            for i in range(3):
                if self.fixed[i]:
                    forces[index,i] = 0.0

    def copy(self):
        return copy.deepcopy(self)


class FixPlane:
    """A constraint class for fixing atomic positions in a plane
    """
    def __init__(self, indices, normal):
        self.indices = indices
        self.normal = np.array(normal) / math.sqrt(np.dot(normal,normal))

    def adjust_positions(self, oldpositions, newpositions):
        for index in self.indices:
            newpositions[index] -= self.normal * \
                np.dot(newpositions[index] - oldpositions[index], self.normal)

    def adjust_forces(self, positions, forces):
        for index in self.indices:
            forces[index] -= self.normal * np.dot(forces[index], self.normal) 

    def copy(self):
        return copy.deepcopy(self)


class FixLine:
    """A constraint class for fixing atomic positions on a line
    """
    def __init__(self, indices, direction):
        self.indices = indices
        self.direction = np.array(direction) / math.sqrt(np.dot(direction,direction))

    def adjust_positions(self, oldpositions, newpositions):
        for index in self.indices:
            tangent = np.dot(newpositions[index] - oldpositions[index], self.direction)
            newpositions[index] = oldpositions[index] + tangent * self.direction

    def adjust_forces(self, positions, forces):
        for index in self.indices:
            tangent = np.dot(forces[index], self.direction)
            forces[index] = tangent * self.direction

    def copy(self):
        return copy.deepcopy(self)

    
class FrictionSimulation:
    """A class defining a simulation environment.

    The class contains convenient methods for setting up a simulation.
    """
    
    # Default values
    timestep = 2.0
    cell_width = 10.0

    def __init__(self):
        self.dynamics = None
        self.system = None
        self.trajectory = None
        self.calc = None
        self.old_positions = None
        self.prev_positions = None
        self.prev_forces = None
        self.work = np.zeros(0,dtype=float)
        self.timestep = FrictionSimulation.timestep
        self.cell_width = FrictionSimulation.cell_width
        self.system = ase.Atoms()
        self.calc = pysic.Pysic()
        self.system.set_calculator(self.calc)
        pysic.utility.error.set_warning_level(1)


    def __del__(self):
        if not self.trajectory is None:
            self.trajectory.close()



    def create_slab(self, element, bottom_z=0.0, top_z=None, z_cells=3, xy_cells=4):
        """Creates a slab of atoms in the FCC (face-centered cubic) lattice structure.
        
        The slab is 2D infinite in xy-plane, its surfaces exposed in the z-direction.
        The structure is FCC, see e.g. http://en.wikipedia.org/wiki/Face_Centered_Cubic
                
        Parameters:
        
        element: string
            the chemical symbol for the atoms of the slab, e.g., 'Au'
        bottom_z: real number
            the z-coordinate of the lower surface of the slab
        top_z: real number
            the z-coordinate of the top surface of the slab
        xy_cells: integer
            the number of lattice cells to be created in x and y directions - as the
            width of the simulation box is determined by the :class:`FrictionSimulation`
            ('cell_width'), this will determine the lattice constant
        z_cells: integer
            the number of lattice cells to be created in z direction
        """

        lattice_constant = self.cell_width / xy_cells
        print "Lattice constant of the created slab: {}".format(lattice_constant)
        slab = bulk(element, 'fcc', a=lattice_constant, cubic=True)*(xy_cells,xy_cells,z_cells)

        if top_z is None:
            slab.translate([0.0, 0.0, bottom_z])
        elif bottom_z == 0.0:
            slab.translate([0.0, 0.0, top_z-lattice_constant*(z_cells-0.5)])
        else:
            print "Only give either the top or the bottom of the slab as an argument."
            print "Setting the z-position accirding to bottom_z = {bz}".format(bz=str(bottom_z))
            slab.translate([0.0, 0.0, bottom_z])

        self.system += slab
        self.system.set_cell([[lattice_constant*xy_cells,0.0,0.0],
                        [0.0,lattice_constant*xy_cells,0.0],
                        [0.0,0.0,100.0]])
        self.system.set_pbc([True, True, False])
        
    def create_random_atoms(self,number,element,bottom_z,top_z, minimum_distance=2.0):
        """Adds the given number of atoms randomly in a given volume.
        
        The atoms are added randomly and for each added atom it is checked that
        it is not too close to previously added atoms. If you try to insert too many
        atoms, there will not be enough room and the algorithm stops after failing
        for long enough.
        
        
        Parameters:
        
        number: integer
            the number of atoms to add
        element: string
            the chemical symbol of the atoms, e.g., 'C'
        bottom_z: real number
            the minimum z coordinate the atoms may get
        top_z: real number
            the maximum z coordinate the atoms may get
        minimum_distance: real number
            the minimum allowed separation between the added atoms
        """
        
        if top_z < bottom_z:
            print "Maximum z must be greater than minimum z. Swapping the values."
            (top_z, bottom_z) = (bottom_z, top_z)            
        
        imps = 0
        tries = 0
        positions = []
        minimum_squared = minimum_distance*minimum_distance
        while imps < number:
            if tries >= 1000:
                print "Cannot find space for new atoms, terminating."
                break
                
            rx = nr.random() * self.cell_width
            ry = nr.random() * self.cell_width
            rz = nr.random() * (top_z-bottom_z) + bottom_z
            new_position = [rx, ry, rz]
        
            if imps == 0:
                # add the first atom
                positions = np.array([new_position])
                imps += 1
            else:
                # check that the new position is not too close to previous ones
                for old_pos in positions:
                
                    separation = np.array([0.0, 0.0, 0.0])
                    # take periodic boundary conditions in to account
                    for i in range(3):
                        separation[i] = new_position[i]-old_pos[i]
                        if separation[i] > self.cell_width/2.0:
                            separation[i] - self.cell_width
                        elif separation[i] < self.cell_width/2.0:
                            separation[i] + self.cell_width

                    distance_squared = np.dot(separation, separation)
                    if distance_squared < minimum_squared:
                        tries += 1
                        break
                else:
                    # if all the previous positions were far away
                    positions = np.append(positions, [new_position], axis=0)
                    imps += 1
                    tries = 0
                    
        print "Creating {ni} atoms.".format(ni=str(imps))
        atoms = Atoms(element+str(imps), positions)
        
        self.system += atoms
        
        self.system.set_cell([ [self.cell_width, 0.0, 0.0],
                         [0.0, self.cell_width, 0.0],
                         [0.0, 0.0, 100.0] ])
        self.system.set_pbc([True, True, False])
        
    def create_atoms(self,element,positions):
        """Creates atoms of the given element in the given positions.

        The atoms are stored in :attr:`tools.FrictionSimulation.pieces`. 
        Once you've added all the pieces you need, finalize the system with
        :meth:`tools.FrictionSimulation.finalize_system`.
        
        Parameters:
        
        element: string
            the chemical symbol of the atoms, e.g., 'C'
        positions: array of real numbers
            the xyz coordinates of the atoms, in an array, e.g., ``[[0, 0, 0], [1, 1, 1]]``
        """
        atoms = Atoms(element+str(len(positions)), np.array(positions))
        
        self.system += atoms
        
        
    def create_interaction(self,element_pair, strength, equilibrium_distance, type='LJ'):
        """Creates an interaction between the two given types of atoms.
        
        The interactions are described by the simple Lennard-Jones potential
        (http://en.wikipedia.org/wiki/Lennard-Jones_potential)
        
        .. math::

           V(r) = \\varepsilon \\left[ \\left( \\frac{\\sigma}{r} \\right)^12 - \\left( \\frac{\\sigma}{r} \\right)^6 \\right]
           
        This is not a realistic interaction except for some very special cases
        such as noble gases. However, it is simple and will work fine as a placeholder.

        The method creates an interaction between two types of atoms as defined by
        ``element_pair``. The types may be the same, so for instance if you have ``'Au'``
        atoms in your system, you can make them interact with 
        ``create_interaction( ['Au', 'Au'], ...)``
        
        Parameters:
        
        element_pair: string list
            the types that will interact, e.g, ``['Au', 'Au']`` or ``['B', 'C']``
        strength: real number
            the depth of the potential energy minimum (:math:`-\\min_r V(r)`) for the pair interaction, in eV
        equilibrium_distance: real number
            the atomic distance, :math:`r`, at which :math:`V(r)` has its minimum
        """
        epsilon = 4*strength
        sigma = equilibrium_distance * 0.890898718
        
        potential = pysic.Potential('LJ', symbols=element_pair, 
                                    cutoff=1.9*sigma, cutoff_margin=0.8*sigma,
                                    parameters=[epsilon, sigma])
      
        self.calc.add_potential(potential)
        
        
    def attach_with_springs(self, index_list1, index_list2, strength, cutoff=5.0):
        """Attaches the atoms from one list with the atoms in another list with springs.
        
        The two lists of atoms need to have an equal number of atoms, :math:`N`, and the 
        method will create harmonic interactions between the listed atoms such that the 
        atoms are connected in order: for ``index_list1 = [0, 1, 2]``, 
        ``index_list2 = [3, 4, 5]`` (:math:`N=3`), 
        the interactions will be created for pairs 0-3, 1-4, and 2-5.
        
        The potential energy of the interaction is
        
        .. math::
        
          V(r) = \\frac{1}{2} k r^2,
          
        where :math:`r` is the distance between the atoms. That means the equilibrium
        distance is 0. You can use the method, for instance, to create harmonic potential
        wells for atoms by attaching them to constrained ghost atoms.
        
        Parameters:
        
        index_list1: integer list
            indices of atoms to be connected
        index_list2: integer list
            indices of atoms list1 is being connected to
        strength: real number
            strength of the attachment, spring constant :math:`k = \\text{strength}/N`
        cutoff: real number
            cutoff for the interaction - if the atoms are farther away from each other than this, the potential is ignored.
            cutoff is restricted to cutoff < cell_width (default 10.0) because otherwise springs would be attached to periodic replicas as well.
        """
        n = len(index_list1)
        if(len(index_list2) != n):
            raise Error("The index lists need to have the same number of atoms.")
        
        k = strength/n
        
        if cutoff >= self.cell_width:
            raise Error("cutoff must be smaller than cell_width")
        
        for atom1, atom2 in zip(index_list1, index_list2):
            newpot = pysic.Potential('spring',indices=[atom1,atom2],
                                     cutoff=cutoff,parameters=[strength, 0.0])
            self.calc.add_potential(newpot)
         
    def continue_from_trajectory(self, filename='simulation.traj', frame=-1):
        """Read the system geometry from a previously saved trajectory.
        
        Parameters:
        
        filename: string
            name of the trajectory file to be read
        frame: integer
            the number of the configuration to be read. by default this is -1, meaning 
            the *last* configuration of the trajectory is read.
        """
    
        traj = PickleTrajectory(filename)
        self.system = traj[frame]
        self.system.set_calculator(self.calc)
        self.remove_constraints()
        
                        
    def list_atoms(self):
        """Lists the elements, indices, and positions of all atoms in the system.
        """
        i = 0        
        for at in self.system:
                print "index = ",i,", element = ",at.symbol,", position = ",at.position
                i += 1
                
    def duplicate_atoms(self, element, indices):
        """Creates a set of atoms with the same positions as the atoms with the given indices.
        
        The new atoms will be of the given element. This function can be used to create
        ghost atoms to be used with :meth:`~friction_tools.FrictionSimulation.attach_with_springs`
        or duplicate parts of the system to be moved in the correct place with
        :meth:`~friction_tools.FrictionSimulation.move_atoms`
        
        Parameters:
        
        element: string
            the chemical symbol of the atoms, e.g., 'C' 
        indices: integer list
            indices of the atoms to be copied
        """
        i = 0   
        newpos = []
        for at in self.system:               
                if i in indices:
                    newpos.append(at.position)
                i += 1
        if len(newpos) > 0:
            self.create_atoms(element=element,positions=np.array(newpos))
        else:
            print "no atoms to duplicate"
                
    def remove_atoms(self, indices):
        """Removes atoms with the given indices.

        Parameters:
        
        indices: integer list
            indices of the atoms to be removed
        """
        for i in range(len(self.system),-1,-1):
            if i in indices:
                del self.system[i]
                
                
    def set_velocities(self, velocity, indices):
        """Gives the atoms initial velocities.
        
        You must give the indices of the atoms whose velocities you are defining as a list.
        If all the atoms should receive the same velocity, ``velocity`` can be a single vector.
        You can also, instead, give a list of velocities to define the velocities one by one.
        In that case, the length of the velocity and index lists must match.
        
        Parameters:
        
        velocity: real vector or a list of vectors
            the velocity or velocities of th atoms
        indices: integer list
            the indices of the atoms affected
        """
        
        if len(velocity) == 3:
            velocity = np.array(len(indices)*[velocity])
        
        if len(velocity) != len(indices):
            print "the velocity and index lists have different lengths in 'set_velocities'"
        
        index = 0
        for i in indices:
            self.system[i].momentum = np.array(velocity[index])*self.system[i].mass
            index += 1
        
    def set_temperature(self, temperature=50):
        """Gives the atoms initial velocities according to the `Maxwell-Boltzmann distribution <http://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution>`_.
        
       The default temperature is 50 K. (Room temperature is about 300 K.)

        Parameters:
        
        temperature: float
            the initial temperature in K - 50 K by default
        """
        if self.system is None:
            print "You need to 'finalize_system' before setting the temperature."
            return
        
        MaxwellBoltzmannDistribution(self.system, temperature*units.kB)
        print ("Set initial velocities according to %6.1f K." % temperature)


    def create_dynamics(self, dt=timestep,
                        temperature=None, coupled_indices=None, strength=0.01):
        """Defines the type of molecular dynamics to be run.

        By default, energy conserving dynamics are run, i.e,
        the microcanonical ensemble is used. That is, the normal
        `Velocity Verlet <http://en.wikipedia.org/wiki/Velocity_Verlet#Velocity_Verlet>`_ 
        algorithm is applied.

        NVT dynamics are run if the argument ``temperature`` is given.         
        (NVT: number of particles, volume, and temperature are kept constant - temperature only approximately.)
        This is done using `Langevin dynamics <http://en.wikipedia.org/wiki/Langevin_dynamics>`_.

        This Langevin thermostat adds a fictious friction force and a random 
        noise force on all atoms. This may lead to computational artifacts if there is collective
        motion. It's especially important not to apply the thermostat on atoms which are being
        dragged or move at a fixed velocity, as the thermostat will add artificial friction
        on them. In order to select the atoms which are being affected by the thermostat,
        you can specify their indices as a list with the ``coupled_indices`` argument.

        - **Note** The atomic structure cannot be altered after this method is called (if temperature is set). Create all the atoms before calling this.

        Parameters:
        
        dt: real number
            timestep :math:`\\Delta t` used in integrating the equations of motion, in fs - 2 fs by default.
        temperature: real number
            temperature in K, if you wish to apply a thermostat
        coupled_indices: integer list
            indices of the atoms the thermostat affects
        strength: real number
            strength of the thermostat, larger number meaning a stronger deviation from Newtonian mechanics - 0.01 by default
        """

        if not temperature is None:
            temp = temperature
                
            if coupled_indices is None:
                lan = strength
            else:       
                lan = []
                for at in range(len(self.system)):
                    if at in coupled_indices:
                        lan.append(strength)
                    else:
                        lan.append(0.0)
                lan = np.array(lan)

            self.dynamics = LangevinMod(self.system, dt*units.fs, temp*units.kB, lan, fixcm=False)
            if pysic.get_cpu_id() == 0:
                  print ("Set up Langevin dynamics")

        else:
                # Create constant energy dynamics
                self.dynamics = VelocityVerlet(self.system, dt*units.fs)
                if pysic.get_cpu_id() == 0:
                    print ("Set up Verlet dynamics")

        self.timestep = dt




    def calculate_potential_energy(self):
        """Returns the potential energy of the system.
        """
        return self.system.get_potential_energy()
        

    def run_simulation(self, time=10, steps=None):
        """Runs molecular dynamics (MD).

        This method executes the actual simulation according to the parameters defined
        earlier by other methods.

        The simulation may take a long time to run, so it is a very good idea to first run a short
        test run to get an idea on the execution time. You can time the simulation with::

          >>> from tools import *
          >>> s = FrictionSimulation()
          >>> # define the system...
          >>> import time
          >>> t0 = time.time()
          >>> s.run_simulation(steps=5)
          >>> t1 = time.time()
          >>> print "The simulation took "+str(t1-t0)+" s."

        The simulation runs for a given length of time, which can be defined as either simulated time or
        number of molecular dynamics steps. If the argument ``steps`` is given, the simulation runs for this
        many steps. If ``time`` is given, the simulation is run for this long (in fs). Even if time :math:`t` is
        specified, the simulation length is still determined in steps :math:`n_\\mathrm{steps} = t/\\Delta t`,
        where :math:`\\Delta t` is the timestep.

        Parameters:

        time: float
            simulation time in fs, 10 fs by default
        steps: integer
            simulation time in MD steps
        """

        if isinstance(self.dynamics,BFGSLineSearch):
            self.dynamics.run(fmax=0.05)
            return

        if steps is None:
            run_steps = int(time/self.timestep)
        else:
            run_steps = steps

        print "Starting a dynamic simulation."

        self.dynamics.run(run_steps)


    def remove_constraints(self):
        """Removes any constraints previously applied on the system.
        """
        self.system.constraints = []


    def fix_positions(self, indices, xyz=[True,True,True]):
        """Constrains the positions of particles.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method freezes coordinates of the atoms whose indices are given.
        By default, the atoms are not allowed to move at all, but the optional argument ``xyz``
        allows one to freeze only some coordinates. For example::
          
          >>> from friction_tools import *
          >>> s = FrictionSimulation()
          >>> # define the system...
          >>> s.fix_positions(indices=[...], xyz=[True,True,False])

        freezes the x and y coordinates of the specified atoms (``xyz = [True,True,False]``).

        Parameters:

        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0
        xyz: list of three booleans
            logical tags specifying if the xyz coordinates are fixed (``True`` means the coordinate is fixed)
        """

        c = FixPositions(indices, xyz)
        self.system.constraints.append(c)
        

    def fix_velocities(self, indices, velocity, xyz=[True,True,True]):
        """Constrains the velocities of particles.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method freezes velocities of the atoms whose indices are given.
        By default, the atoms move with constant velocity, but the optional argument ``xyz``
        allows one to fix only a given component of force. For example::
          
          >>> from tools import *
          >>> s = FrictionSimulation()
          >>> # define the system...
          >>> s.fix_velocities(indices=[...], velocity=[0.01,0.0,0.0], xyz=[True,False,True])

        fixes the velocities in x and z directions to :math:`v_x = 0.01, v_z = 0.0` Å/fs
        (``xyz = [True,False,True]``) but allows the velocity to change freely in the y direction.

        Parameters:

        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0
        velocity: list of three floats
            the velocity given to the particles
        xyz: list of three booleans
            logical tags specifying if the xyz coordinates are fixed (``True`` means the coordinate is fixed)
        """
        c = FixVelocities(indices, velocity, self.timestep, xyz)
        self.system.constraints.append(c)

        
    def fix_positions_on_line(self, indices, direction):
        """Constrains the positions of particles on a line.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method constrains coordinates of the atoms whose indices are given on a line.
        The line is defined by the given direction :math:`\\mathbf{u}` and the initial position
        of the particle :math:`\\mathbf{r}_0`.

        .. math::

          \\mathbf{r} - \\mathbf{r}_0 \\parallel \\mathbf{u}

        Parameters:

        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0
        direction: list of three floats
            the direction ``[ux, uy, uz]`` in which the particles are allowed to move
        """
        c = FixLine(indices, direction)
        self.system.constraints.append(c)


    def fix_positions_on_plane(self, indices, normal):
        """Constrains the positions of particles on a plane.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method constrains coordinates of the atoms whose indices are given on a plane.
        The plane is defined by the given normal :math:`\\mathbf{n}` and the initial position
        of the particle :math:`\\mathbf{r}_0`.

        .. math::

          \\mathbf{r} - \\mathbf{r}_0 \\perp \\mathbf{n}

        Parameters:

        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0
        normal: list of three floats
            the normal vector ``[nx, ny, nz]`` to the plane in which the particles are allowed to move
        """
        c = FixPlane(indices, normal)
        self.system.constraints.append(c)


    def add_constant_force(self, indices, force):
        """Adds a constant external force on atoms.

        Adds an external force on the specified atoms
 
        .. math::

          \\mathbf{F} = f_x \\mathbf{i} + f_y \\mathbf{j} + f_z \\mathbf{k}.

        Parameters:

        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0
        force: list of three floats
            the force ``[fx, fy, fz]`` to be applied on the given atoms

        """
        index_list = []
        for index in indices:
            index_list.append([index])
        pot = pysic.Potential('force', indices=index_list, parameters=force)
        self.calc.add_potential(pot)


    def print_stats(self, time_it=False):
        """Prints the energy and temperature information in human readable format.

        Parameters:

        time_it: boolean
            if ``True``, the calculation is timed and the execution time is printed
        """
        t0 = time.time()
        epot = self.system.get_potential_energy() / len(self.system)
        t1 = time.time()
        ekin = self.system.get_kinetic_energy() / len(self.system)
        t2 = time.time()
        stress = self.system.get_stress()
        t3 = time.time()
        step = self._get_step()

        if pysic.get_cpu_id() == 0:
            if step > 0:
                print ("Timestep: %7i" % step)
            print ("Energy per atom: Epot = %.3f eV,  Ekin = %.3f eV,  Etot = %.3f eV" %
                   (epot, ekin, epot+ekin))
            print ("Temperature = %.3f K\n" %
                   (ekin/(1.5*units.kB)))
            if time_it:
                print ("Energy calculation took %.3f s, force and stress calculation took %.3f s\n" %
                       (t1-t0, t3-t2))


    def write_positions_to_file(self, filename='positions.txt', indices=None, addStepNumber=False):
        """Writes the coordinates of the atoms in a file.

        The format is::

           (atom 0 x) (atom 0 y) (atom 0 z)
           (atom 1 x) (atom 1 y) (atom 1 z)
           (atom 2 x) (atom 2 y) (atom 2 z)
           ...

        
        Parameters:

        filename: string
            name of the file where the information is written
        indices: list of integers
           the indices of the atoms whose info is to be written - note that Python indexing starts from 0. By default all atoms are included
        addStepNumber: boolean
           meant for internal use - if ``True``, adds the MD step number in the filename to prevent overwriting
        """
        if indices is None:
            indices = range(len(self.system))

        output_lines = ""
        for index in indices:
            [x,y,z] = self.system[index].get_position()
            output_lines += (" %20.5f  %20.5f  %20.5f \n" % (x,y,z) )

        filing = add_number_to_filename(filename, self._get_step(), addStepNumber)
        write_file(filing, output_lines)


    def write_velocities_to_file(self, filename='velocities.txt', indices=None, addStepNumber=False):
        """Writes the instantaneous velocities of the atoms in a file.

        - **Note** Due to a bug/feature in ASE, if you fix the velocities of atoms with :meth:`tools.FrictionSimulation.fix_velocities`, the velocities for these atoms are reported to be 0.0. The method :meth:`tools.FrictionSimulation.gather_average_velocities_during_simulation` does report the correct average velocities for all atoms though.

        The format is::

           (atom 0 v_x) (atom 0 v_y) (atom 0 v_z)
           (atom 1 v_x) (atom 1 v_y) (atom 1 v_z)
           (atom 2 v_x) (atom 2 v_y) (atom 2 v_z)
           ...
        
        Parameters:

        filename: string
            name of the file where the information is written
        indices: list of integers
           the indices of the atoms whose info is to be written - note that Python indexing starts from 0. By default all atoms are included
        addStepNumber: boolean
           meant for internal use - if ``True``, adds the MD step number in the filename to prevent overwriting
        """

        if indices is None:
            indices = range(len(self.system))

        output_lines = ""
        for index in indices:
            try:
                [x,y,z] = self.system[index].get_momentum()/self.system[index].get_mass()
            except:
                [x,y,z] = [0.0,0.0,0.0]
            output_lines += (" %20.5f  %20.5f  %20.5f \n" % (x,y,z) )

        filing = add_number_to_filename(filename, self._get_step(), addStepNumber)
        write_file(filing, output_lines)

    def get_avg_velocity(self, indices):
        """ Average velocity of indices
        """

        sum_velocity = [0, 0, 0]
        for index in indices:
            try:
                veloc = (
                    self.system[index].momentum /
                    self.system[index].mass
                )
            except:
                veloc = [0.0, 0.0, 0.0]
            sum_velocity = [sum_velocity[i] + veloc[i] for i in range(0, 3)]
        return [sum_velocity[i] / len(indices) for i in range(0, 3)]

    def write_forces_to_file(self, filename='forces.txt', indices=None, addStepNumber=False):
        """Writes the instantaneous forces of the atoms in a file.

        The format is::

           (atom 0 f_x) (atom 0 f_y) (atom 0 f_z)
           (atom 1 f_x) (atom 1 f_y) (atom 1 f_z)
           (atom 2 f_x) (atom 2 f_y) (atom 2 f_z)
           ...

        
        Parameters:

        filename: string
            name of the file where the information is written
        indices: list of integers
           the indices of the atoms whose info is to be written - note that Python indexing starts from 0. By default all atoms are included
        addStepNumber: boolean
           meant for internal use - if ``True``, adds the MD step number in the filename to prevent overwriting
        """
        
        if indices is None:
            indices = range(len(self.system))

        output_lines = ""
        forces = self.system.get_forces(apply_constraint=False)
        for index in indices:
            [x,y,z] = forces[index]
            output_lines += (" %20.5f  %20.5f  %20.5f \n" % (x,y,z) ) 

        filing = add_number_to_filename(filename, self._get_step(), addStepNumber)
        write_file(filing, output_lines)


    def write_average_force_to_file(self, filename='force_sum.txt', indices=None):
        """Writes the average of instantaneous forces of the atoms in a file.

        The format is::

           (f_x) (f_y) (f_z)
           ...
        
        Parameters:

        filename: string
            name of the file where the information is written
        indices: list of integers
           the indices of the atoms whose info is to be written - note that Python indexing starts from 0. By default all atoms are included
        """
        
        if indices is None:
            indices = range(len(self.system))

        output_lines = ""
        forces = self.system.get_forces(apply_constraint=False)
        force_sum = np.array([0.0,0.0,0.0])
        for index in indices:
            force_sum += forces[index]

        [x,y,z] = force_sum / float(len(indices))
        output_line = (" %20.5f  %20.5f  %20.5f \n" % (x,y,z) )
        append_file(filename, output_line)




    def write_average_position_to_file(self, filename='avr_position.txt', indices=None):
        """Writes the mean of the coordinates of the atoms in a file.

        The format is::

           (x) (y) (z)
        
        Parameters:

        filename: string
            name of the file where the information is written
        indices: list of integers
           the indices of the atoms whose info is to be written - note that Python indexing starts from 0. By default all atoms are included
        """
        
        if indices is None:
            indices = range(len(self.system))

        output_lines = ""
        pos = self.system.get_positions()
        pos_sum = np.array([0.0,0.0,0.0])
        for index in indices:
            pos_sum += pos[index]

        [x,y,z] = pos_sum / float(len(indices))
        output_line = (" %20.5f  %20.5f  %20.5f \n" % (x,y,z) )
        append_file(filename, output_line)





    def write_energy_and_temperature_to_file(self, filename='energies.txt'):
        """Appends the instantaneous energy and temperature information in a file.

        The format is::

           (potential energy) (kinetic energy) (total energy) (temperature)

        Since this method does not write atom-by-atom information, it does not create a new file but
        appends its output at the end of the specified file.

        
        Parameters:

        filename: string
            name of the file where the information is written
        """

        e_pot = self.system.get_potential_energy()
        e_kin = self.system.get_kinetic_energy()
        t = e_kin/(len(self.system)*(1.5*units.kB))
        output_lines = (" %20.5f  %20.5f  %20.5f  %20.5f \n" % (e_pot, e_kin, (e_pot+e_kin), t) )
        append_file(filename, output_lines)


    def _record_positions(self):
        """For internal use. Saves the positions of the atoms."""
        self.old_positions = self.system.get_positions()


    def _record_work(self):
        """For internal use. Integrates work done on all atoms."""
        forces = self.system.get_forces(apply_constraint=False)            
        move = self.system.get_positions()
        move -= self.prev_positions

        for index in range(len(self.system)):
            self.work[index] += np.dot((forces[index]+self.prev_forces[index])*0.5, move[index])
        self.prev_positions = self.system.get_positions()        
        self.prev_forces = forces


    def _clear_work(self):
        """For internal use. Clears the recorded work integrals."""
        self.work = np.zeros(len(self.system),dtype=float)
        self.prev_positions = self.system.get_positions()
        self.prev_forces = self.system.get_forces(apply_constraint=False)


    def write_work_to_file(self, filename='work.txt', indices=None, addStepNumber=False):
        """For internal use. Writes the work done on the atoms in a file.

        The format is::

           (atom 0 work) 
           (atom 1 work) 
           (atom 2 work) 
           ...

        
        Parameters:

        filename: string
            name of the file where the information is written
        indices: list of integers
           the indices of the atoms whose info is to be written - note that Python indexing starts from 0. By default all atoms are included
        addStepNumber: boolean
           meant for internal use - if ``True``, adds the MD step number in the filename to prevent overwriting
        """
        
        if indices is None:
            indices = range(len(self.system))

        output_lines = ""
        
        for index in indices:
            output_lines += (" %20.5f \n" % (self.work[index]) )

        #self._clear_work()
        filing = add_number_to_filename(filename, self._get_step(), addStepNumber)
        write_file(filing, output_lines)


    def write_total_work_to_file(self,filename='total_work.txt', indices=None):
        """For internal use. Writes the work done on the atoms in a file.
        
        Parameters:

        filename: string
            name of the file where the information is written
        indices: list of integers
           the indices of the atoms whose info is to be written - note that Python indexing starts from 0. By default all atoms are included
        """
                
        if indices is None:
            indices = range(len(self.system))

        output_lines = ""
        work_sum = 0.0
        for index in indices:
            work_sum += self.work[index]

        output_line = (" %20.5f \n" % (work_sum) )
        append_file(filename, output_line)


    def write_average_velocities_to_file(self, dt, filename='avr_velocities.txt',
                                         indices=None, addStepNumber=False):
        """For internal use. Writes the average velocities of the atoms in a file.

        The format is::

           (atom 0 v_x) (atom 0 v_y) (atom 0 v_z)
           (atom 1 v_x) (atom 1 v_y) (atom 1 v_z)
           (atom 2 v_x) (atom 2 v_y) (atom 2 v_z)
           ...

        
        Parameters:

        filename: string
            name of the file where the information is written
        indices: list of integers
           the indices of the atoms whose info is to be written - note that Python indexing starts from 0. By default all atoms are included
        addStepNumber: boolean
           meant for internal use - if ``True``, adds the MD step number in the filename to prevent overwriting
        """
        
        if indices is None:
            indices = range(len(self.system))

        output_lines = ""
        for index in indices:
            pos = self.system[index].get_position()
            pos -= self.old_positions[index]
            pos /= dt*units.fs
            [x,y,z] = pos
            
            output_lines += (" %20.5f  %20.5f  %20.5f \n" % (x,y,z) )

        self._record_positions()
        filing = add_number_to_filename(filename, self._get_step(), addStepNumber)
        append_file(filing, output_lines)

        
    def print_stats_during_simulation(self, interval=10):
        """Makes the simulator print energies and temperature during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.print_stats` periodically
        during the molecular dynamics run, every ``interval`` steps. This is a convenient way to
        monitor the progress of your simulation.
        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        """

        self.dynamics.attach(self.print_stats, interval)

    def gather_positions_during_simulation(self, interval=10, filename='positions.txt',
                                           indices=None):
        """Makes the simulator write coordinates of atoms during simulation.

        The format is::

           (atom 0 x) (atom 0 y) (atom 0 z)
           (atom 1 x) (atom 1 y) (atom 1 z)
           (atom 2 x) (atom 2 y) (atom 2 z)
           ...
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_positions_to_file` periodically
        during the molecular dynamics run, every ``interval`` steps.

        The given filename is automatically modified to include the current MD step number.
        For instance, if the filename is given as ``data.txt``, the data from, say,
        step 20 will be written to file ``data00020.txt``.

        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0. By default, all atoms are included.
        """

        self.dynamics.attach(self.write_positions_to_file, interval, filename, indices, True)


    def gather_velocities_during_simulation(self, interval=10, filename='velocities.txt',
                                            indices=None):
        """Makes the simulator write instantaneous velocities of atoms during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`
        - **Note** Due to a bug/feature in ASE, if you fix the velocities of atoms with :meth:`tools.FrictionSimulation.fix_velocities`, the velocities for these atoms are reported to be 0.0. The method :meth:`tools.FrictionSimulation.gather_average_velocities_during_simulation` does report the correct average velocities for all atoms though.

        The format is::

           (atom 0 v_x) (atom 0 v_y) (atom 0 v_z)
           (atom 1 v_x) (atom 1 v_y) (atom 1 v_z)
           (atom 2 v_x) (atom 2 v_y) (atom 2 v_z)
           ...

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_velocities_to_file`
        periodically during the molecular dynamics run, every ``interval`` steps.

        The given filename is automatically modified to include the current MD step number.
        For instance, if the filename is given as ``data.txt``, the data from, say,
        step 20 will be written to file ``data00020.txt``.

        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0. By default, all atoms are included.
        """

        self.dynamics.attach(self.write_velocities_to_file, interval, filename, indices, True)


    def gather_average_velocities_during_simulation(self, interval=10, filename='avr_velocities.txt',
                                                    indices=None):
        """Makes the simulator write average velocities of atoms during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        The format is::

           (atom 0 v_x) (atom 0 v_y) (atom 0 v_z)
           (atom 1 v_x) (atom 1 v_y) (atom 1 v_z)
           (atom 2 v_x) (atom 2 v_y) (atom 2 v_z)
           ...

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_average_velocities_to_file`
        periodically during the molecular dynamics run, every ``interval`` steps.

        The average velocities
        recorded are the averages since the previous recording.
        So, when writing at intervals of :math:`n_\\mathrm{avr}`,
        the average velocity of atom :math:`i` at time :math:`t` will be the atomic displacement
        :math:`\\Delta \\mathbf{r}_i(t) = \\mathbf{r}_i(t) - \\mathbf{r}_{i}(t-\\Delta t)`
        divided by the elapsed time :math:`\\Delta t = n_\\mathrm{avr} \\Delta t'`, with :math:`\\Delta t'` being
        the MD timestep:

        .. math::
           
              \\mathbf{v}_{i,\\mathrm{avr}}(t-\\Delta t, t) = \\frac{\\Delta \\mathrm{r}_i}{\\Delta t_i}.
           
        The given filename is automatically modified to include the current MD step number.
        For instance, if the filename is given as ``data.txt``, the data from, say,
        step 20 will be written to file ``data00020.txt``.
        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0. By default, all atoms are included.
        """
        self._record_positions()
        self.dynamics.attach(self.write_average_velocities_to_file, interval,
                             interval*self.timestep, filename, indices, True)

    def gather_forces_during_simulation(self, interval=100, filename='forces.txt',
                                        indices=None):
        """Makes the simulator write instantaneous forces of atoms during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        The format is::

           (atom 0 f_x) (atom 0 f_y) (atom 0 f_z)
           (atom 1 f_x) (atom 1 f_y) (atom 1 f_z)
           (atom 2 f_x) (atom 2 f_y) (atom 2 f_z)
           ...

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_forces_to_file`
        periodically during the molecular dynamics run, every ``interval`` steps.

        The given filename is automatically modified to include the current MD step number.
        For instance, if the filename is given as ``data.txt``, the data from, say,
        step 20 will be written to file ``data00020.txt``.

        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0. By default, all atoms are included.
        """

        self.dynamics.attach(self.write_forces_to_file, interval, filename, indices, True)


    def gather_average_position_during_simulation(self, interval=100, filename='avr_position.txt',
                                            indices=None):
        """Makes the simulator write averages of positions of atoms during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        The format is::

           (x) (y) (z)
           

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_average_position_to_file`
        periodically during the molecular dynamics run, every ``interval`` steps.

        The data is written in a single file where each line represents a point in the time series.
        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0. By default, all atoms are included.
        """

        write_file(filename, "") # clear the file
        self.dynamics.attach(self.write_average_position_to_file, interval, filename, indices)



    def gather_average_force_during_simulation(self, interval=100, filename='avr_force.txt',
                                            indices=None):
        """Makes the simulator write averages of instantaneous forces of atoms during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        The format is::

           (f_x) (f_y) (f_z)
           

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_average_force_to_file`
        periodically during the molecular dynamics run, every ``interval`` steps.

        The data is written in a single file where each line represents a point in the time series.
        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0. By default, all atoms are included.
        """

        write_file(filename, "") # clear the file
        self.dynamics.attach(self.write_average_force_to_file, interval, filename, indices)


    def gather_total_work_during_simulation(self, interval=100, filename='total_work.txt',
                                            indices=None):
        """Makes the simulator write the total work done on atoms during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_total_work_to_file`
        periodically during the molecular dynamics run, every ``interval`` steps.

        The data is written in a single file where each line represents a point in the time series.
        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0. By default, all atoms are included.
        """

        write_file(filename, "") # clear the file
        self._clear_work()
        self.dynamics.attach(self._record_work, interval=1)
        self.dynamics.attach(self.write_total_work_to_file, interval, filename, indices)



    def gather_energy_and_temperature_during_simulation(self, interval=100, filename='energy.txt'):
        """Makes the simulator write energy and temperature during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        The format is::

          (potential energy)  (kinetic energy)  (total energy)  (temperature)

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_energy_and_temperature_to_file`
        periodically during the molecular dynamics run, every ``interval`` steps.

        The data is written in a single file where each line represents a point in the time series.
        
        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        """
        write_file(filename, "") # clear the file
        self.dynamics.attach(self.write_energy_and_temperature_to_file, interval, filename)
        

    def gather_work_during_simulation(self, interval=100, filename='work.txt', indices=None):
        """Makes the simulator write work done on atoms during simulation.
        
        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        The format is::

          (atom 0 work)
          (atom 1 work)
          (atom 2 work)
          ...

        This method instructs the simulator to call :meth:`tools.FrictionSimulation.write_work_to_file`
        periodically during the molecular dynamics run, every ``interval`` steps.

        Work is defined as the line integral of the force. The recorded values are the work done
        on each individual atom since the previous recording.

        So, when writing at intervals of :math:`n_\\mathrm{avr}`, the work on atom :math:`i`
        at time :math:`t` is

        .. math::

           W_i(t-\\Delta t, t) = \\int_{t-\\Delta t}^t \\mathbf{F}_i \\cdot \\mathrm{d}\\mathbf{r}_i(t)

        with :math:`\\Delta t = n_\\mathrm{avr} \\Delta t'`, :math:`\\Delta t'` being
        the MD timestep.
        
        The given filename is automatically modified to include the current MD step number.
        For instance, if the filename is given as ``data.txt``, the data from, say,
        step 20 will be written to file ``data00020.txt``.

        
        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        indices: list of integers
            the indices of the atoms to be constrained - note that Python indexing starts from 0. By default, all atoms are included.
        """
        self._clear_work()
        self.dynamics.attach(self._record_work, interval=1)
        self.dynamics.attach(self.write_work_to_file, interval, filename, indices, True)
        

    def save_trajectory_during_simulation(self, interval=100, filename='simulation.traj'):
        """Saves the simulation trajectory in binary format.

        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method saves the simulation as a trajectory file, from which the positions of
        particles can be extracted at the different points of time of the simulation.

        
        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written
        """
        
        if not self.trajectory is None:
            #self.trajectory.close()
            mode = 'a'
        else:
            mode = 'w'

        self.trajectory = PickleTrajectory(filename, mode, self.system, master=( pysic.get_cpu_id() == 0 ))
        self.dynamics.attach(self.trajectory.write, interval)


    def take_snapshots_during_simulation(self, interval=100, filename='snapshot.png'):
        """Takes png-format snapshots of the system during the simulation.

        - **Note** call this method before running the simulation with :meth:`tools.FrictionSimulation.run_simulation`

        This method generates images of the system in png format.

        Parameters:

        interval: integer
            the number of MD steps between consecutive information writing
        filename: string
            name of the file where the data is written        
        """
        self.dynamics.attach(self.take_snapshot, interval, filename, True)

    
    def load_state_from_trajectory(self, step=-1, filename='simulation.traj'):
        """Loads the geometry from a saved trajectory.

        If you have a saved trajectory from :meth:`tools.FrictionSimulation.save_trajectory_during_simulation`,
        this method can be used for reading in a state from the trajectory.

        
        Parameters:

        step: integer
            the timestep which is restored - by default the last recorded state is chosen
        filename: string
            name of the trajectory file
        """
        traj_in = PickleTrajectory(filename)
        self.system = traj_in[step]

    
    def _get_step(self):
        """Returns the number of steps the dynamics has run.
        """
        try:
            return self.dynamics.get_number_of_steps()
        except:
            return 0



    def take_snapshot(self, filename='snapshot.png', addStepNumber=False):
        """Makes a png format snapshot of the system.

        Parameters:

        filename: string
            name of the trajectory file
        addStepNumber: boolean
           meant for internal use - if ``True``, adds the MD step number in the filename to prevent overwriting
        """
        
        filing = add_number_to_filename(filename, self._get_step(), addStepNumber)
        ase.io.write(filing, self.system, format='png', rotation='-80x')


    def write_xyzfile(self, filename='system.xyz'):
        """Makes an xyz format geometry file of the system.

        XYZ is a commonly used, very simple geometry file format which includes
        information on the types of atoms and their coordinates. This method writes
        such a file based on the current simulation geometry.

        XYZ files can be read with almost any atomic visualization tool.
        
        Parameters:

        filename: string
            name of the trajectory file
        """
        ase.io.write(filename, self.system, format='xyz')

    def get_indices_z_less_than(self, z_limit):
        """Returns the indices of the atoms whose z-coordinate is less than the given value, as a list.

        This can be used for easily finding, e.g., the bottom of a slab.

        
        Parameters:

        z_limit: float
            the limiting z value
        """
        return [atom.index for atom in self.system if atom.z < z_limit]

    def get_indices_z_more_than(self, z_limit):
        """Returns the indices of the atoms whose z-coordinate is more than the given value, as a list.

        This can be used for easily finding, e.g., the top of a slab.

        
        Parameters:

        z_limit: float
            the limiting z value
        """
        return [atom.index for atom in self.system if atom.z > z_limit]

    def get_indices_by_element(self, element):
        """Returns the indices of atoms of given element.
        
        Parameters:
                
        element: string
            the chemical symbol of the atoms, e.g., 'C'
        """
        return [atom.index for atom in self.system if atom.symbol == element]

def read_arrays_from_files(filename, starting_index, ending_index, index_step):
    """Reads numeric arrays from files.

    If you have run a simulation and produced files containing numeric data through the
    methods ``gather_xxx_during_simulation``, this function can be used for reading the
    data in Python.

    You should have files with names such as ``positions00100.txt`` etc. To read them in,
    give the base name ``positions.txt`` as the ``filename``, and the range of integers through the
    other arguments: ``starting_index`` and ``ending_index`` are the first and last index to be read,
    respectively, while ``index_step`` is the increment.

    Example::

      >>> read_arrays_from_files('data.txt',10,100,10)

    will read the files::

      data00010.txt
      data00020.txt
      ...
      data00100.txt

    The function returns a list containing all the arrays extracted from the files.

    
    Parameters:

    filename: string
        the base filename to be read
    starting_index: integer
        the first index to be read
    ending_index: integer
        the last index to be read
    index_step: integer
        the increment in indices
    """
    array_collection = []
    for index in range(starting_index, ending_index+index_step, index_step):
        filing = add_number_to_filename(filename, index, True)
        new_array = read_array_from_file(filing)
        array_collection.append(new_array)

    return array_collection

def read_array_from_file(filename):
    """Parses a file and returns a numeric array of the data found.

    Parameters:

    filename: string
        name of the file to be read
    """
    array = np.fromfile(filename, dtype=float, count=-1, sep=' ')
    f=file(filename)
    lines = f.readlines()
    f.close()
    array.shape = (len(lines), array.size/len(lines))
    return array

def write_array_to_file(filename, array):
    """Writes an array to a file.

    **Note** The method calls the Numpy routine::

      array.tofile(filename, sep=" ", format="%20.5f")

    which results in all the values stored in the array being
    written in one line. If you read it back in with
    :meth:`tools.read_array_from_file`, you get a 1D array
    containing the values. To fix the shape, specify it explicitly::

     >>> arry = [[1, 2], [3, 4]]
     array([[1, 2],
            [3, 4]])
     >>> write_array_to_file(filename='array.txt', arry)
     >>> new_arry = read_array_from_file('array.txt')
     >>> new_arry
     array([[1, 2, 3, 4]])
     >>> new_arry.shape = (2,2)
     >>> new_arry
     array([[1, 2],
            [3, 4]])

    Parameters:

    filename: string
        name of the file to be written
    array: numeric array
        the array to be written
    """
    array.tofile(filename, sep=" ", format="%20.5f")

def add_arrays_by_component(array1, array2):
    """Adds two arrays component by component.

    Equal to ``np.array(array1) + np.array(array2)``

    Parameters:

    array1: numeric array of floats
         an array
    array2: numeric array of floats
         an array
    """
    return np.array(array1) + np.array(array2)

def subtract_arrays_by_component(array1, array2):
    """Subtracts two arrays component by component.

    Equal to ``np.array(array1) - np.array(array2)``

    Parameters:

    array1: numeric array of floats
         an array
    array2: numeric array of floats
         an array
    """
    return np.array(array1) - np.array(array2)

def multiply_arrays_by_component(array1, array2):
    """Multiplies two arrays component by component.

    Equal to ``np.array(array1) * np.array(array2)``

    Parameters:

    array1: numeric array of floats
         an array
    array2: numeric array of floats
         an array
    """
    return np.array(array1) * np.array(array2)

def divide_arrays_by_component(array1, array2):
    """Divides two arrays component by component.

    Equal to ``np.array(array1) / np.array(array2)``

    Parameters:

    array1: numeric array of floats
         an array
    array2: numeric array of floats
         an array
    """
    return np.array(array1) / np.array(array2)

def add_to_array(array, real):
    """Adds a real number to an array component by component.

    Equal to ``np.array(array) + real``

    Parameters:

    array: numeric array of floats
         an array
    real: float
         a number
    """
    return np.array(array) + real

def multiply_array(array, real):
    """Multiplies an array with a real number component by component.

    Equal to ``np.array(array) * real``

    Parameters:

    array: numeric array of floats
         an array
    real: float
         a number
    """
    return np.array(array1) * real

def create_array(dimensions, value=0):
    """Creates an array of the given size.

    Parameters:

    dimensions: list of integers
        dimensions of the array
    value: float
        the initial value for the elements of the array
    """
    return np.zeros(dimensions) + value

def concatenate_arrays_by_rows(array1, array2):
    """Joins two arrays together by appending row by row.

    Parameters:

    array1: numeric array of floats
         an array
    array2: numeric array of floats
         an array
    """
    cols = array1.shape[1]
    rows1 = array1.shape[0]
    rows2 = array2.shape[0]
    newarray = np.zeros((rows1+rows2,cols))
    for i in range(rows1):
        newarray[i,:] = array1[i,:]
    for i in range(rows2):
        newarray[i+rows1,:] = array2[i,:]
    return newarray

def concatenate_arrays_by_columns(array1, array2):
    """Joins two arrays together by appending column by column.

    Parameters:

    array1: numeric array of floats
         an array
    array2: numeric array of floats
         an array
    """
    cols1 = array1.shape[1]
    rows = array1.shape[0]
    cols2 = array2.shape[1]
    newarray = np.zeros((rows,cols1+cols2))
    for i in range(cols1):
        newarray[:,i] = array1[:,i]
    for i in range(cols2):
        newarray[:,i+cols1] = array2[:,i]
    return newarray


def shift_array_by_rows(array, steps=1, wrap=False):
    """Moves the elements in an array by rows.

    Parameters:

    array: numeric array of floats
         an array
    steps: integer
         number of rows the elements are moved
    wrap: boolean
         if ``True``, elements that go over the bounds of the array appear on the other side
    """
    newarray = np.zeros_like(array)
    newarray[steps:] = np.array(array)[:-steps]
    if wrap:
        newarray[:steps] = np.array(array)[-steps:]
    return newarray
    
def shift_array_by_columns(array, steps=1, wrap=False):
    """Moves the elements in an array in by columns.

    Parameters:

    array: numeric array of floats
         an array
    steps: integer
         number of rows the elements are moved
    wrap: boolean
         if ``True``, elements that go over the bounds of the array appear on the other side
    """
    newarray = np.zeros_like(array)
    newarray[:,steps:] = np.array(array)[:,:-steps]
    if wrap:
        newarray[:,:steps] = np.array(array)[:,-steps:]
    return newarray


def crop_array(top, left, bottom, right, array):
    """Crops an array (i.e., returns a subarray).

    The indices for cropping should be given as arguments.
    The order is ``top, left, bottom, right``.
    For ``top`` and ``left``, they are the first indices to be included,
    while ``bottom`` and ``right`` are the first indices to be dropped.

    For instance, if you have the array::

      arry = [[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9]]

    The top left :math:`2 \\times 2` array is obtained with::

      >>> crop_array(0, 0, 2, 2, arry)
      [[1, 2],
       [4, 5]]


    The command is equivalent to::

      >>> array.copy()[top:bottom,left:right]

    Parameters:

    top: integer
        first index from the top to be included
    left: integer
        first integer from the left to be included
    bottom: integer
        first index from the top to be left out
    rigth: integer
        first index from the left to be left out
    """
    return array.copy()[top:bottom,left:right]


def calculate_row_sums(array):
    """Sums all the rows of a 2D array and returns a row 1D array containing the results.
    
    Parameters:

    array: numeric array
        the array to be manipulated
    """
    if len(array.shape) < 3:
        height = array.shape[0]
        
        result = np.zeros(height)
    
        for i in range(height):
            row = array[i,:]
            row_sum = np.sum(row)
            result[i] = row_sum

    else:
        raise Exception("Row sum is only defined for 1- or 2-dimensional arrays")

    return result

def calculate_row_norms(array):
    """Calculates euclidian norms for all the rows of a 2D array and returns a row 1D array containing the results.
    
    Parameters:

    array: numeric array
        the array to be manipulated
    """
    if len(array.shape) < 3:
        height = array.shape[0]
        
        result = np.zeros(height)
    
        for i in range(height):
            row = array[i,:]
            row_sum = sqrt(np.dot(row,row))
            result[i] = row_sum

    else:
        raise Exception("Row norm is only defined for 1- or 2-dimensional arrays")

    return result

def calculate_row_maxima(array):
    """Finds the maxima for all the rows of a 2D array and returns a row 1D array containing the results.
    
    Parameters:

    array: numeric array
        the array to be manipulated
    """
    if len(array.shape) < 3:
        height = array.shape[0]
        
        result = np.zeros(height)
    
        for i in range(height):
            row = array[i,:]
            row_sum = np.max(row)
            result[i] = row_sum

    else:
        raise Exception("Row max is only defined for 1- or 2-dimensional arrays")

    return result

def calculate_row_minima(array):
    """Finds the minima for all the rows of a 2D array and returns a row 1D array containing the results.
    
    Parameters:

    array: numeric array
        the array to be manipulated
    """
    if len(array.shape) < 3:
        height = array.shape[0]
        
        result = np.zeros(height)
    
        for i in range(height):
            row = array[i,:]
            row_sum = np.min(row)
            result[i] = row_sum

    else:
        raise Exception("Row min is only defined for 1- or 2-dimensional arrays")

    return result

def calculate_column_sums(array):
    """Sums all the columns of a 2D array and returns a row 1D array containing the results.
    
    Parameters:

    array: numeric array
        the array to be manipulated
    """
    if len(array.shape) == 2:
        width = array.shape[1]
        
        result = np.zeros(width)
    
        for i in range(width):
            column = array[:,i]
            column_sum = np.sum(column)
            result[i] = column_sum
    elif len(array.shape) == 1:
        result = array

    else:
        raise Exception("Column sum is only defined for 1- or 2-dimensional arrays")

    return result


def calculate_column_norms(array):
    """Calculates euclidian norms for all the columns of a 2D array and returns a row 1D array containing the results.
    
    Parameters:

    array: numeric array
        the array to be manipulated
    """
    if len(array.shape) == 2:
        width = array.shape[1]
        
        result = np.zeros(width)
    
        for i in range(width):
            column = array[:,i]
            column_sum = sqrt(np.dot(column,column))
            result[i] = column_sum
    elif len(array.shape) == 1:
        result = array

    else:
        raise Exception("Column norm is only defined for 1- or 2-dimensional arrays")

    return result

def calculate_column_maxima(array):
    """Finds the maxima for all the columns of a 2D array and returns a row 1D array containing the results.
    
    Parameters:

    array: numeric array
        the array to be manipulated
    """
    if len(array.shape) == 2:
        width = array.shape[1]
        
        result = np.zeros(width)
    
        for i in range(width):
            column = array[:,i]
            column_sum = np.max(column)
            result[i] = column_sum
    elif len(array.shape) == 1:
        result = array

    else:
        raise Exception("Column max is only defined for 1- or 2-dimensional arrays")

    return result

def calculate_column_minima(array):
    """Finds the minima for all the columns of a 2D array and returns a row 1D array containing the results.
    
    Parameters:

    array: numeric array
        the array to be manipulated
    """
    if len(array.shape) == 2:
        width = array.shape[1]
        
        result = np.zeros(width)
    
        for i in range(width):
            column = array[:,i]
            column_sum = np.min(column)
            result[i] = column_sum
    elif len(array.shape) == 1:
        result = array

    else:
        raise Exception("Column min is only defined for 1- or 2-dimensional arrays")

    return result


def plot_array(x_values, y_values=None, y_error=None, x_range=None, y_range=None):
    """Plots data.

    If only one column of data is given, it is plotted against a running index.

    If two columns of data are given, the first is treated as x and the second as y.

    **Note** You should also use :func:`tools.hold_plots`.

    Parameters:

    x_values: list of floats
        If it's the only input given, it is treated as the y-variable. If ``y_values`` is given, this is the x-variable.
    y_values: list of floats
        Values of the y-variable.
    y_error: list of floats
        Values for errorbars of y.
    x_range: list of two floats
        Plotting range for x
    y_range: list of two floats
        Plotting range for y
    """
    if y_values == None:
        xs = np.array(range(len(x_values)))
        ys = np.array(x_values)
    else:
        xs = np.array(x_values)
        ys = np.array(y_values)
    if not x_range is None:        
        plt.ylim(x_range[0],x_range[1])
    if not y_range is None:        
        plt.ylim(y_range[0],y_range[1])
    if y_error == None:
        plt.plot(xs,ys)
    else:
        es = np.array(y_error)
        errorbar(xs, ys, yerr=es, fmt='ro')

def plot_zero_line():
    """Adds a horizontal line at y=0 to the plot.
    """
    plt.axhline()

def hold_plots():
    """Hold the plot visible.

    By default, your plots will vanish from the screen after plotting. This
    function prevents it. The point is, you can call :func:`tools.plot_array` several
    times and the generated plots will be added together and shown when you call
    this function.
    """
    plt.show()

def trajectory_to_xyz(filename_in='simulation.traj',filename_out='simulation.xyz'):
    """Converts a trajectory file to an XYZ file.

    Parameters:

    filaname_in: string
        name of the trajectory file
    filaname_out: string
        name of the XYZ file
    """
    # Write the trajectory to an xyz file
    if pysic.get_cpu_id() == 0:
        traj = PickleTrajectory(filename_in)
        movie = []
        for snapshot in traj:
            movie.append(snapshot)
        ase.io.write(filename_out,movie,format='xyz')

def add_number_to_filename(filename, number, act=True):
    if act:
        name = filename.split('.')
        try:
            filing = ( (name[0]+"%5.5i."+name[1]) % number )
        except:
            filing = ( (filename+"%5.5i") % number )
    else:
        filing = filename

    return filing

def write_file(filename, output_lines):
    """Writes the given string in a file.

    Note that you should provide only a single string. If you want to split the output
    to several lines, use the newline character ``\\n``.

    If the file exists, it is overwritten.

    Parameters:

    filename: string
        name of the file to be written
    output_lines: string
        the string to be written
    """
    if pysic.get_cpu_id() == 0: # if we run MPI, only write with the master CPU
        f = file(filename,'w')
        f.write(output_lines)
        f.close()

def append_file(filename, output_lines):
    """Writes the given string in a file.

    Note that you should provide only a single string. If you want to split the output
    to several lines, use the newline character '\\n'.

    This function does not overwrite a file, but appends the output at the end.

    Parameters:

    filename: string
        name of the file to be written
    output_lines: string
        the string to be written
    """

    if pysic.get_cpu_id() == 0: # if we run MPI, only write with the master CPU
        f = file(filename,'a')
        f.write(output_lines)
        f.close()


