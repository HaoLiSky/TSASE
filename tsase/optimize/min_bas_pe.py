import numpy as np
import os
from ase import io, units
from ase.optimize.optimize import Dynamics
from ase.optimize.fire import FIRE
from ase.optimize import QuasiNewton
from tsase.optimize.sdlbfgs import SDLBFGS
from ase.units import kB
from ase.parallel import world
from ase.io.trajectory import PickleTrajectory
from ase.md import VelocityVerlet
from ase.md import MDLogger
import tsase
import sys

class Hopping(Dynamics):
    """Basin hopping algorithm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116

    and 

    David J. Wales and Harold A. Scheraga, Science, Vol. 285, 1368 (1999)
    """

    def __init__(self, atoms,
		 temperature=100, # K, initial temperature 
                 optimizer=SDLBFGS,
                 fmax=0.1,
                 dr=0.4,
                 logfile='-', 
                 trajectory=None,
                 optimizer_logfile='-',
                 local_minima_trajectory='local_minima.con',
                 adjust_cm=True,
                 mss=0.05,
                 minenergy=None,
                 distribution='uniform',
                 adjust_step_size=None,
                 adjust_every = None,
                 target_ratio = 0.5,
                 adjust_fraction = 0.05,
                 significant_structure = True,  # displace from minimum at each move
                 significant_structure2 = False, # displace from global minimum found so far at each move
                 pushapart = 0.4,
                 jumpmax=None,
                 jmp = 7, # number of jump steps taken in BHOJ
    		 molecular_dynamics = False, # trial move using md when true
		 dimer_a = 0.001,
		 dimer_d = 0.01,
		 dimer_steps = 20,
		 timestep = 0.1, # moleculare dynamics time step
		 mdmin = 2, # number of minima to pass in md before stopping
		 history_weight = 0.0, # the weight factor of history >= 0 
                 adjust_temp = False, # dynamically adjust the temperature in BH acceptance
		 minimaHopping_acceptance = False, # use MH acceptance criteria instead of BH
                 minimaHopping_history = True, # use history in MH acceptance criteria
		 beta1 = 1.04, # temperature adjustment parameter
		 beta2 = 1.04, # temperature adjustment parameter
		 beta3 = 1.0 / 1.04, # temperature adjustment parameter
	         Ediff0 = 0.5,  # eV, initial energy acceptance threshold
        	 alpha1 = 0.98,  # energy threshold adjustment parameter
        	 alpha2 = 1. / 0.98,  # energy threshold adjustment parameter
		 minima_threshold = 2  # round energies to how many decimal places
                 ):
        Dynamics.__init__(self, atoms, logfile, trajectory)
	self.temperature = temperature
        self.optimizer = optimizer
        self.fmax = fmax
        self.dr = dr
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
        self.adjust_step = adjust_step_size
        self.adjust_every = adjust_every
        self.target_ratio = target_ratio
        self.adjust_fraction = adjust_fraction
        self.significant_structure = significant_structure
        self.significant_structure2 = significant_structure2
        self.pushapart = pushapart
        self.jumpmax = jumpmax
        self.jmp = jmp
	self.dimer_a = dimer_a
	self.molecular_dynamics = molecular_dynamics
	self.dimer_d = dimer_d
	self.dimer_steps = dimer_steps
	self.timestep = timestep
	self.mdmin = mdmin
	self.w = history_weight
        self.num_accepted_steps = 0
        self.adjust_temp = adjust_temp
	self.mh_accept = minimaHopping_acceptance
	self.mh_history = minimaHopping_history
	self.beta1 = beta1
	self.beta2 = beta2
	self.beta3 = beta3
	self.Ediff = Ediff0
	self.alpha1 = alpha1
	self.alpha2 = alpha2
	self.minima_threshold = minima_threshold

        # when a MD sim. has passed a local minimum:
        self.passedminimum = PassedMinimum()

        self.mss = mss
	# dictionary for found local minima
        # keys will the potential energy rounded to self.minima_threshold digits left of the decimal
        # values will be number of times the potential energy has been visited 
	self.minima = {}
        self.initialize()

    def initialize(self):
        self.positions = 0.0 * self.atoms.get_positions()
        self.Emin = self.get_energy(self.atoms.get_positions()) or 1.e32 
        self.rmin = self.atoms.get_positions()
        self.positions = self.atoms.get_positions()
        self.call_observers()
        self.log(-1, self.Emin, self.Emin,self.dr)
                
    def log(self, step, En, Emin,dr):
        if self.logfile is None:
            return
        name = self.__class__.__name__
        self.logfile.write('%s: step %d, energy %15.6f, emin %15.6f'
                           % (name, step, En, Emin))
	if not self.molecular_dynamics:
            self.logfile.write(', dr %15.6f'
                               % (dr))
	self.logfile.write(', Temperature %12.4f' % (self.temperature))
        if self.mh_accept:
            self.logfile.write(', Ediff %12.4f\n'
                           % (self.Ediff))
        else:
            self.logfile.write(', w %12.4f\n'
       	                   % (self.w))
        self.logfile.flush()

    def find_ene_match(self, En, Eo):
        """ determines if En is the same PE as any previously visited minima."""
        approxEn = round(En,self.minima_threshold)
        approxEo = round(Eo,self.minima_threshold)
        try:
            countEo = self.minima[approxEo]
        except KeyError:
            # approxEo isn't in the dict of pervious local min
            # this should only be a possiblity at step 0 of the run
            countEo = 0
        # check if approxEn is a previously visited minimum
        if approxEn in self.minima:
            countEn = self.minima[approxEn]
            if approxEn == approxEo:
                return 1, countEn, countEo
            return 2, countEn, countEo
        # approxEn is a new local minimum
        else:
            return None, 0, countEo

    def update_minima(self, En, Eo):
        """Update the dictionary of local minima to include this new location
           and return True if En is a new local minima."""
        self.num_accepted_steps += 1
        approxEn = round(En,self.minima_threshold)
        approxEo = round(Eo,self.minima_threshold)
        if approxEn in self.minima:
            self.minima[approxEn] += 1
            return False, approxEn, approxEo
        else:
            self.minima[approxEn] = 1
            return True, approxEn, approxEo

    def _maxwellboltzmanndistribution(self,masses,N,temp,communicator=world):
        # For parallel GPAW simulations, the random velocities should be
        # distributed.  Uses gpaw world communicator as default, but allow
        # option of specifying other communicator (for ensemble runs)
        xi = np.random.standard_normal((len(masses), 3))
        if N is not None:
            xi = N
        momenta = xi * np.sqrt(masses * temp)[:, np.newaxis]
        communicator.broadcast(xi, 0)
        return momenta

    def MaxwellBoltzmannDistribution(self,N,temp,communicator=world,
                                     force_temp=False):
        """Sets the momenta to a Maxwell-Boltzmann distribution. temp should be
        fed in energy units; i.e., for 300 K use temp=300.*units.kB. If
        force_temp is set to True, it scales the random momenta such that the
        temperature request is precise.
        """
	momenta = self._maxwellboltzmanndistribution(self.atoms.get_masses(),N,temp,
                                                communicator)
        self.atoms.set_momenta(momenta)
        if force_temp:
            temp0 = self.atoms.get_kinetic_energy() / len(self.atoms) / 1.5
            gamma = temp / temp0
            self.atoms.set_momenta(self.atoms.get_momenta() * np.sqrt(gamma))

    def _molecular_dynamics(self, step, N):
        """Performs a molecular dynamics simulation, until mdmin is
        exceeded. If resuming, the file number (md%05i) is expected."""
        mincount = 0
        energies, oldpositions = [], []
        thermalized = False
        if not thermalized:
            self.MaxwellBoltzmannDistribution(N,
                                         temp=self.temperature * kB,
                                         force_temp=True)
        if (step > 1):
            os.remove('md.log')
            os.remove('md.traj')
        traj = io.Trajectory('md.traj', 'a', self.atoms)
        dyn = VelocityVerlet(self.atoms, dt=self.timestep * units.fs)
        log = MDLogger(dyn, self.atoms, 'md.log',
                       header=True, stress=False, peratom=False)
        dyn.attach(log, interval=1)
        dyn.attach(traj, interval=1)
        while mincount < self.mdmin:
            dyn.run(1)
            energies.append(self.atoms.get_potential_energy())
            passedmin = self.passedminimum(energies)
            if passedmin:
                mincount += 1
            oldpositions.append(self.atoms.positions.copy())
        # Reset atoms to minimum point.
        self.atoms.positions = oldpositions[passedmin[0]]

    def get_minimum(self):
        """Return minimal energy and configuration."""
        self.atoms.set_positions(self.rmin)
        print 'get_minimum',self.Emin
        return self.Emin

    def get_energy(self, positions):
        """Return the energy of the nearest local minimum."""
        if np.sometrue(self.positions != positions):
            self.positions = positions
            self.atoms.set_positions(positions)
 
            try:
                if self.optimizer == QuasiNewton:
                    # the mss we use for SDLBFGS makes QN go crazy :(
		    opt = self.optimizer(self.atoms,
				         logfile=self.optimizer_logfile)
		else:
                    opt = self.optimizer(self.atoms,
                                         logfile=self.optimizer_logfile)
                #                        maxstep=self.mss)
                #    opt = self.optimizer(self.atoms, 
                #                     logfile=self.optimizer_logfile,
                #                     maxstep=self.mss)
                opt.run(fmax=self.fmax)
		# this caused an error when using QN local optimizer because
		# energy is not an attriute of the class Hopping
                #self.energy = self.atoms.get_potential_energy()
                self.local_min_pos = self.atoms.get_positions()
            except:
                # Something went wrong.
                # In GPAW the atoms are probably to near to each other.
                return None
        
        return self.atoms.get_potential_energy()

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

    def move(self, step, ro):
        """Move atoms by a random step."""
        if not self.molecular_dynamics:
            if self.distribution == 'uniform':
                disp = np.random.uniform(-self.dr, self.dr, (len(self.atoms), 3))
            elif self.distribution == 'gaussian':
                disp = np.random.normal(0,self.dr,size=(len(self.atoms), 3))
            elif self.distribution == 'linear':
                distgeo = self.get_dist_geo_center()
                disp = np.zeros(np.shape(self.atoms.get_positions()))
                for i in range(len(disp)):
                    maxdist = self.dr*distgeo[i]
                    disp[i] = np.random.uniform(-maxdist,maxdist,3)
            elif self.distribution == 'quadratic':
                distgeo = self.get_dist_geo_center()
                disp = np.zeros(np.shape(self.atoms.get_positions()))
                for i in range(len(disp)):
                    maxdist = self.dr*distgeo[i]*distgeo[i]
                    disp[i] = np.random.uniform(-maxdist,maxdist,3)
            else:
                disp = np.random.uniform(-1*self.dr, self.dr, (len(self.atoms), 3))
            if self.significant_structure == True:
                rn = self.local_min_pos + disp
            elif self.significant_structure2 == True:
                ro  = self.get_minimum()
                rn = ro + disp
            else:
                rn = ro + disp
            rn = self.push_apart(rn)
         #   atoms.set_positions(rn)
            self.atoms.set_positions(rn)
            if self.cm is not None:
                cm = self.atoms.get_center_of_mass()
                self.atoms.translate(self.cm - cm)
        else :
            # Move atoms using Dimer and an MD step
            dimer = ModifiedDimer()
            N = dimer(self.atoms, self.dimer_a, self.dimer_d, self.dimer_steps)
            self._molecular_dynamics(step, N)
        rn = self.atoms.get_positions()
        world.broadcast(rn, 0)
        self.atoms.set_positions(rn)
        return self.atoms.get_positions()

    def acceptance_MH(self, steps, ro, Eo, maxtemp):
        """Adjusts parameters and positions based  on minima hopping acceptance criteria."""
        rejectnum = 0
        for step in range(steps):
            positionsOld = self.atoms.get_positions()
            En = None
            self.steps += 1
            rn = self.move(step,ro)
            En = self.get_energy(rn)
            while En is None:
                rn = self.move(step,ro)
                En = self.get_energy(rn)
            if En < self.Emin:
                self.Emin = En
                self.rmin = self.atoms.get_positions()
                self.call_observers()
            self.log(step, En, self.Emin,self.dr)
            match, approxEn, approxEo = self.find_ene_match(En, Eo)
            if match is not None:
                if match == 1:
                # re-found last minimum
                    self.temperature *= self.beta1
                elif self.mh_history:
                # re-found previously found minimum
                    self.temperature *= self.beta2
            else:
            # must have found a new minimum
                self.temperature *= self.beta3
            accept = False
            if (En < (Eo + self.Ediff)):
                self.Ediff *= self.alpha1
                accept = True
            else:
                self.Ediff *= self.alpha2
                rejectnum += 1
            if rejectnum > self.jumpmax:
                for i in range(0,self.jmp):
                    rn = self.move(step,rn)
                accept = True
            if accept:
                rejectnum = 0
                ro = rn
                Eo = En
                self.update_minima(En, Eo)
                if self.lm_trajectory is not None:
                    tsase.io.write_con(self.lm_trajectory,self.atoms,w='a')
            else:
                self.atoms.set_positions(positionsOld)
            if self.minenergy is not None:
                if Eo < self.minenergy:
                    break
            if maxtemp and maxtemp < self.temperature:
                  break

    def acceptance_BH(self, steps, ro, Eo):
        """Boltzmann acceptance criteria for basin hopping."""
        acceptnum = 0
        recentaccept = 0
        rejectnum = 0
        for step in range(steps):
            positionsOld = self.atoms.get_positions()
            En = None
            rn = None
            self.steps += 1
            while En is None:
                rn = self.move(step,ro)
                En = self.get_energy(rn)
            if En < self.Emin:
                self.Emin = En
                self.rmin = self.atoms.get_positions()
                self.call_observers()
            self.log(step, En, self.Emin,self.dr)
            match, countEn, countEo = self.find_ene_match(En, Eo)
            if match is not None and self.adjust_temp:
                if match == 1:
                # re-found last minimum
                    self.temperature *= self.beta1
                else:
                # re-found previously found minimum
                    self.temperature *= self.beta2
            elif self.adjust_temp:
            # must have found a new minimum
                self.temperature *= self.beta3
            accept = False;
            if Eo >= En and self.w == 0.0:
                accept = True
            else:
                totalMin = len(self.minima)
                hn = 0
                ho = 0
                if self.num_accepted_steps > 0:
                    if match is not None:
                            hn = countEn / float(self.num_accepted_steps)
                    ho = countEo / float(self.num_accepted_steps)
                kT = self.temperature * kB
                # occasionally overflowing the floating point :(
                #accept = np.exp(((Eo - En) + (self.w * (ho - hn))) / kT) > np.random.uniform()
                val = ((Eo - En) + (self.w * (ho - hn))) / kT
                if val < 1.0:
                    accept = np.exp(val) > np.random.uniform()
                else: # accept the new position
                    accept = True;
            if rejectnum > self.jumpmax:
                #JMP???
                for i in range(0,self.jmp):
                    rn = self.move(step,rn)
                accept = True
            if accept:
                acceptnum += 1.
                recentaccept += 1.
                rejectnum = 0
                if self.significant_structure2 == True:
                    ro = self.local_min_pos.copy()
                else:
                    ro = rn.copy()
                Eo = En
                if self.lm_trajectory is not None:
                    tsase.io.write_con(self.lm_trajectory,self.atoms,w='a')
            else:
                rejectnum += 1
                self.atoms.set_positions(positionsOld)
            if self.minenergy != None:
                if Eo < self.minenergy:
                    break
            if self.adjust_step == True:
                if step % self.adjust_every == 0:
                    ratio = float(acceptnum)/float(self.steps)
                    ratio = float(recentaccept)/float(self.adjust_every)
                    recentaccept = 0.
                    if ratio > self.target_ratio:
                       self.dr = self.dr * (1+self.adjust_fraction)
                    elif ratio < self.target_ratio:
                        self.dr = self.dr * (1-self.adjust_fraction)

    def run(self, steps, maxtemp = None):
        """Hop the basins for defined number of steps."""
        self.steps = 0
        ro = self.positions
        Eo = self.get_energy(ro)
        if self.mh_accept:
            self.acceptance_MH(steps, ro, Eo, maxtemp)
        else:
            self.acceptance_BH(steps, ro, Eo)
        self.get_minimum()

    def get_dist_geo_center(self):
        position = self.atoms.get_positions()
        geocenter = np.sum(position,axis=0)/float(len(position))
        distance = np.zeros(len(position))
        for i in range(len(distance)):
            vec = position[i]-geocenter
            distance[i] = np.sqrt(np.vdot(vec,vec))
        distance /= np.max(distance)  #np.sqrt(np.vdot(distance,distance))
        return distance


    def _unique_minimum_position(self):
        """Identifies if the current position of the atoms, which should be
        a local minima, has been found before."""
        unique = True
        dmax_closest = 99999.
        compare = CompareEnergies()
        #self._read_minima()
        for minimum in self._minima:
            dmax = compare(minimum, self._atoms)
            if dmax < self._minima_threshold:
                unique = False
            if dmax < dmax_closest:
                dmax_closest = dmax
        return unique, dmax_closest


class ModifiedDimer:
# Class that moves the initial velocity vector of a MD escape trial
# toward a direction with low curvature
    def __call__(self,atoms,dimer_a,dimer_d,dimer_steps):
        p = atoms.copy()
        #count = counter
        oldPos = atoms.get_positions()
        lj = tsase.calculators.lj(cutoff=35.0)
        p.set_calculator(lj)
        N = self.gradientDimer(p,dimer_a,dimer_d,dimer_steps)
        atoms.set_positions(oldPos)
        return N

    def perpForce(self,p, N):
        # function that calculates the perpendicular force
        force = p.get_forces()
        f = force - np.vdot(force,N) * N
        return f

    def escapeDirection(self,x, y):
        # function that calculates the escape direction vector
        diff = y - x
	return diff / np.linalg.norm(diff)

    def gradientDimer(self,p,dimer_a,dimer_d,dimer_steps):
        localOpt = tsase.optimize.SDLBFGS(p,logfile = None)
        localOpt.run()
        x = p.get_positions()
        a = dimer_a
        d = dimer_d
        # random uniform vector
        N = np.random.uniform(-1, 1, (len(p), 3))
        N = N / np.linalg.norm(N)
        y = x + d * N
        p.set_positions(y)
        maxIteration = dimer_steps
        iteration = 0

        # After a few steps the iteration is stopped before a locally
        # optimal lowest curvature mode is found
        while iteration < maxIteration:
            f = self.perpForce(p,N)
            y = y + a * f
            N = self.escapeDirection(x,y)
            y = x + d * N
            p.set_positions(y)
            iteration += 1
        return N

class PassedMinimum:
    """Simple routine to find if a minimum in the potential energy surface
    has been passed. In its default settings, a minimum is found if the
    sequence ends with two downward points followed by two upward points.
    Initialize with n_down and n_up, integer values of the number of up and
    down points. If it has successfully determined it passed a minimum, it
    returns the value (energy) of that minimum and the number of positions
    back it occurred, otherwise returns None."""

    def __init__(self, n_down=2, n_up=2):
        self._ndown = n_down
        self._nup = n_up

    def __call__(self, energies):
        if len(energies) < (self._nup + self._ndown + 1):
            return None
        status = True
        index = -1
        for i_up in range(self._nup):
            if energies[index] < energies[index - 1]:
                status = False
            index -= 1
        for i_down in range(self._ndown):
            if energies[index] > energies[index - 1]:
                status = False
            index -= 1
        if status:
            return (-self._nup - 1), energies[-self._nup - 1]


class CompareEnergies:
    """Class that compares the potential energy of 'M' and 'Mcurrent'
    or 'M' with all other local minima perviously found
    """

    def __call__(self, atoms1, atoms2):
        atoms1 = atoms1.copy()
        atoms2 = atoms2.copy()
        dmax = self._indistinguishable_compare(atoms1, atoms2)
        return dmax

    def _indistinguishable_compare(self, atoms1, atoms2):
        lj = tsase.calculators.lj()
        atoms1.set_calculator(lj)
        atoms2.set_calculator(lj)
        difference = atoms2.get_potential_energy() - atoms1.get_potential_energy()
        dmax = np.absolute(difference)
        return dmax


