import numpy as np

from ase.optimize.optimize import Dynamics
from tsase.optimize.sdlbfgs import SDLBFGS
from ase.units import kB
from ase.parallel import world
from ase.io.trajectory import PickleTrajectory
import tsase
import sys

class BasinHopping(Dynamics):
    """Basin hopping algorithm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116

    and 

    David J. Wales and Harold A. Scheraga, Science, Vol. 285, 1368 (1999)
    """

    def __init__(self, atoms,
                 temperature=100 * kB,
                 optimizer=SDLBFGS,
                 fmax=0.1,
                 dr=0.1,
                 active_ratio = 1.0, #define the number of atoms moved each time
                 logfile='-', 
                 trajectory=None,
                 optimizer_logfile='-',
                 local_minima_trajectory='local_minima.con',
                 adjust_cm=True,
                 mss=0.2,
                 minenergy=None,
                 distribution='uniform',
                 adjust_step_size = None,
                 target_ratio = 0.5,
                 adjust_fraction = 0.05,
                 significant_structure = False,  # displace from minimum if accept
                 significant_structure2 = False, # displace from global minimum found so far at each move
                 pushapart = 0.4,
                 jumpmax=None
                 ):
        Dynamics.__init__(self, atoms, logfile, trajectory)
        self.kT = temperature
        self.optimizer = optimizer
        self.fmax = fmax
        self.dr = dr
        self.active_ratio = active_ratio

        if adjust_cm:
            self.cm = atoms.get_center_of_mass()
        else:
            self.cm = None

        self.optimizer_logfile = optimizer_logfile
        self.lm_trajectory = local_minima_trajectory
        if isinstance(local_minima_trajectory, str):
            tsase.io.write_con(self.lm_trajectory,atoms,w='w')
        self.minenergy = minenergy
        self.distribution = distribution
        self.adjust_step_zie = adjust_step_size
        self.target_ratio = target_ratio
        self.adjust_fraction = adjust_fraction
        self.significant_structure = significant_structure
        self.significant_structure2 = significant_structure2
        self.pushapart = pushapart
        self.jumpmax = jumpmax
        self.mss = mss
        self.initialize()

    def initialize(self):
        self.positions = 0.0 * self.atoms.get_positions()
        self.Emin = self.get_energy(self.atoms.get_positions()) or 1.e32 
        self.rmin = self.atoms.get_positions()
        self.positions = self.atoms.get_positions()
        self.local_min_pos = self.atoms.get_positions()
        self.call_observers()
        self.log(-1, self.Emin, self.Emin,self.dr)
                
    def run(self, steps):
        """Hop the basins for defined number of steps."""
        self.steps = 0
        ro = self.positions
        Eo = self.get_energy(ro)
        acceptnum = 0
        recentaccept = 0
        rejectnum = 0
        for step in range(steps):
            En = None
            self.steps += 1
            while En is None:
                rn = self.move(ro)
                En = self.get_energy(rn)
            if En < self.Emin:
                self.Emin = En
                self.rmin = self.atoms.get_positions()
                self.call_observers()
            self.log(step, En, self.Emin,self.dr)
            if Eo >= En:
                accept = True
            else:
                accept = np.exp((Eo - En) / self.kT) > np.random.uniform()
            if rejectnum > self.jumpmax:
                accept = True
                rejectnum = 0
            if accept:
                acceptnum += 1.
                recentaccept += 1.
                rejectnum = 0
                #Lei: update structure with accepted local minimum
                #if self.significant_structure2 == True:
                #    ro = self.local_min_pos.copy()
                if self.significant_structure == True:
                    ro = self.local_min_pos.copy()
                else:
                    ro = rn.copy()
                Eo = En
                if self.lm_trajectory is not None:
                    tsase.io.write_con(self.lm_trajectory,self.atoms,w='a')
            else:
                rejectnum += 1
            if self.minenergy != None:
                if Eo < self.minenergy:
                    break
            #Lei: merge two parameters 'adjust_every' and 'adjust_step_size' as one
            if self.adjust_step_size in not None:
                if step % self.adjust_step_size == 0:
                    ratio = float(acceptnum)/float(self.steps)
                    ratio = float(recentaccept)/float(self.adjust_step_size)
                    recentaccept = 0.
                    if ratio > self.target_ratio:
                       self.dr = self.dr * (1+self.adjust_fraction)
                    elif ratio < self.target_ratio:
                        self.dr = self.dr * (1-self.adjust_fraction)

    def log(self, step, En, Emin,dr):
        if self.logfile is None:
            return
        name = self.__class__.__name__
        self.logfile.write('%s: step %d, energy %15.6f, emin %15.6f, dr %15.6f\n'
                           % (name, step, En, Emin,dr))
        self.logfile.flush()

    def move(self, ro):
        """Move atoms by a random step."""
        atoms = self.atoms
        disp = np.zeros(np.shape(atoms.get_positions()))
        while np.alltrue(disp == np.zeros(np.shape(atoms.get_positions()))):
            if self.distribution == 'uniform':
                disp = np.random.uniform(-self.dr, self.dr, (len(atoms), 3))
            elif self.distribution == 'gaussian':
                disp = np.random.normal(0,self.dr,size=(len(atoms), 3))
            elif self.distribution == 'linear':
                distgeo = self.get_dist_geo_center()
                disp = np.zeros(np.shape(atoms.get_positions()))
                for i in range(len(disp)):
                    maxdist = self.dr*distgeo[i]
                #    disp[i] = np.random.normal(0,maxdist,3)
                    disp[i] = np.random.uniform(-maxdist,maxdist,3)
            elif self.distribution == 'quadratic':
                distgeo = self.get_dist_geo_center()
                disp = np.zeros(np.shape(atoms.get_positions()))
                for i in range(len(disp)):
                    maxdist = self.dr*distgeo[i]*distgeo[i]
                #    disp[i] = np.random.normal(0,maxdist,3)
                    disp[i] = np.random.uniform(-maxdist,maxdist,3)
            else:
                disp = np.random.uniform(-1*self.dr, self.dr, (len(atoms), 3))

            #Lei: set all other disp to zero except those selected to move
            #     the number of atoms that can be moved is defined by int(active_ratio * len(atoms))
            if self.active_ratio is not None:
               fix_space = len(atoms) - int(self.active_ratio * len(atoms))
               fix_atoms = random.sample(range(len(atoms)), fix_space)
               for i in range(len(fix_atoms)):
                   disp[fix_atoms[i]] = (0.0, 0.0, 0.0)

        #Lei: suppose to only update 'ro' with local minimum when accept == true
        #if self.significant_structure == True:
        #    rn = self.local_min_pos + disp
        if self.significant_structure2 == True:
            ro,reng = self.get_minimum()
            rn = ro + disp
        else:
            rn = ro + disp
        rn = self.push_apart(rn)
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
        print 'get_minimum',self.Emin
        return self.Emin, atoms

    def get_energy(self, positions):
        """Return the energy of the nearest local minimum."""
        if np.sometrue(self.positions != positions):
            self.positions = positions
            self.atoms.set_positions(positions)
 
            try:
               #Lei: enable 'FIRE' optimizer
               if self.optimizer.__name__ == "FIRE":
                  opt = self.optimizer(self.atoms,
                                       maxmove = 1.0,
                                       dt = 0.2, dtmax = 1.0,
                                           logfile=self.optimizer_logfile)
               else:
                  opt = self.optimizer(self.atoms,
                                           logfile=self.optimizer_logfile,
                                           maxstep=self.mss)
                #    opt = self.optimizer(self.atoms, 
                #                     logfile=self.optimizer_logfile,
                #                     maxstep=self.mss)
                opt.run(fmax=self.fmax)
                self.energy = self.atoms.get_potential_energy()
                self.local_min_pos = self.atoms.get_positions()
            except:
                # Something went wrong.
                # In GPAW the atoms are probably to near to each other.
                return None
        
        return self.energy

    def push_apart(self,positions):
        movea = np.zeros(np.shape(positions))
        alpha = 0.025
        for w in range(500):
            moved = 0
            movea = np.zeros(np.shape(positions))
            for i in range(len(positions)):
                for j in range(i+1,len(positions)):
                    d = positions[i] - positions[j]
                    magd = np.sqrt(np.vdot(d,d))
                    if magd < self.pushapart:
                        moved += 1
                        vec = d/magd
                        movea[i] += alpha *vec
                        movea[j] -= alpha *vec
            positions += movea
            if moved == 0:
                break
        return positions

    def get_dist_geo_center(self):
        position = self.atoms.get_positions()
        geocenter = np.sum(position,axis=0)/float(len(position))
        distance = np.zeros(len(position))
        for i in range(len(distance)):
            vec = position[i]-geocenter
            distance[i] = np.sqrt(np.vdot(vec,vec))
        distance /= np.max(distance)  #np.sqrt(np.vdot(distance,distance))
        return distance 
