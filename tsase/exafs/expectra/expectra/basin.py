import numpy as np

from ase.optimize.optimize import Dynamics
from ase.optimize.fire import FIRE
from ase.units import kB
from ase.parallel import world
from ase.io.trajectory import Trajectory


class BasinHopping(Dynamics):
    """Basin hopping algorithm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116

    and 

    David J. Wales and Harold A. Scheraga, Science, Vol. 285, 1368 (1999)
    """

    def __init__(self, atoms,
                 opt_calculator = None,
                 exafs_calculator = None,
                 alpha = 0,
                 temperature=100 * kB,
                 optimizer=FIRE,
                 fmax=0.1,
                 dr=0.1,
                 logfile='-', 
                 trajectory='lowest.traj',
                 optimizer_logfile='-',
                 local_minima_trajectory='local_minima.traj',
                 adjust_cm=True):
        """Parameters:

        atoms: Atoms object
            The Atoms object to operate on.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.
        pseudo_pot: pseudo potential defined
        """
        Dynamics.__init__(self, atoms, logfile, trajectory)
        self.kT = temperature
        self.optimizer = optimizer
        self.fmax = fmax
        self.dr = dr
        self.opt_calculator = opt_calculator
        self.exafs_calculator = exafs_calculator
        self.alpha = alpha

        if adjust_cm:
            self.cm = atoms.get_center_of_mass()
        else:
            self.cm = None

        self.optimizer_logfile = optimizer_logfile
        self.lm_trajectory = local_minima_trajectory
        if isinstance(local_minima_trajectory, str):
            self.lm_trajectory = Trajectory(local_minima_trajectory,
                                                  'w', atoms)

        self.initialize()

    def initialize(self):
        self.positions = 0.0 * self.atoms.get_positions()
        self.Umin = self.get_energy(self.atoms.get_positions()) or 1.e32
        self.chi_deviation = 100
        self.rmin = self.atoms.get_positions()
        self.positions = self.atoms.get_positions()
        self.call_observers()
        self.log(-1, self.Umin, self.chi_deviation, self.Umin,self.Umin)
                
    def run(self, steps):
        """Hop the basins for defined number of steps."""

        ro = self.positions
        Eo = self.get_energy(ro)
        chi_devi_o = self.get_chi_deviation(ro)
        Uo = Eo + alpha * chi_devi_o

        for step in range(steps):
            Un = None
            while Un is None:
                rn = self.move(ro)
                En = self.get_energy(rn)
                chi_devi_n = self.get_chi_deviation(rn)
                Un = En + alpha * chi_devi_n

            if Un < self.Umin:
                # new minimum found
                self.Umin = Un
                self.rmin = self.atoms.get_positions()
                self.call_observers()
            self.log(step, En, chi_devi_n, Un, self.Umin)

            #accept or reject?
            accept = np.exp((Uo - Un) / self.kT) > np.random.uniform()
            if accept:
                ro = rn.copy()
                Uo = Un

    def log(self, step, En, chi_devi_n, Un, Umin):
        if self.logfile is None:
            return
        name = self.__class__.__name__
        self.logfile.write('%s: step %d, energy %15.6f, chi_deviation %15.6f, pseudoPot %15.6f, emin %15.6f\n'
                           % (name, step, En, chi_devi_n, Un, Umin))
        self.logfile.flush()

    def move(self, ro):
        """Move atoms by a random step."""
        atoms = self.atoms
        # displace coordinates
        disp = np.random.uniform(-1., 1., (len(atoms), 3))
        rn = ro + self.dr * disp
        atoms.set_positions(rn)
        if self.cm is not None:
            cm = atoms.get_center_of_mass()
            atoms.translate(self.cm - cm)
        rn = atoms.get_positions()
        world.broadcast(rn, 0)
        atoms.set_positions(rn)
        return atoms.get_positions()

    def get_minimum(self):
        """Return minimal energy and configuration."""
        atoms = self.atoms.copy()
        atoms.set_positions(self.rmin)
        return self.Umin, atoms

    def get_energy(self, positions):
        """Return the energy of the nearest local minimum."""
        if np.sometrue(self.positions != positions):
            self.positions = positions
            self.atoms.set_positions(positions)

            #opt_calculator can be any calculator compatible with ASE
            if self.opt_calculator is not None:
                try:
                    self.atoms.set_calculator(self.opt_calculator)
                    self.energy = self.atoms.get_potential_energy()
                except:
                    # Something went wrong.
                    # In GPAW the atoms are probably to near to each other.
                    return None

                return self.energy

            #run only if opt_calculator is not set
            try:
                opt = self.optimizer(self.atoms, 
                                     logfile=self.optimizer_logfile)
                opt.run(fmax=self.fmax)
                if self.lm_trajectory is not None:
                    self.lm_trajectory.write(self.atoms)

                self.energy = self.atoms.get_potential_energy()
            except:
                # Something went wrong.
                # In GPAW the atoms are probably to near to each other.
                return None

            return self.energy

        return self.energy

    def get_chi_deviation(self, positions):
        """Return the standard deviation of chi between calculated and
        experimental."""
        if np.sometrue(self.positions != positions):
            self.positions = positions
            self.atoms.set_positions(positions)

            try:
                self.atoms.set_calculator(self.exafs_calculator)
                self.chi_deviation = self.atoms.get_potential_energy()
            except:
                # Something went wrong.
                # In GPAW the atoms are probably to near to each other.
                return None
   
        return self.chi_deviation

