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
#import atoms as geometry
import eon.fileio as fileio
#import write_con as write

class Hopping(Dynamics):
    """Basin hopping algorithm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116

    and 

    David J. Wales and Harold A. Scheraga, Science, Vol. 285, 1368 (1999)
    """

    def __init__(self, atoms,
		 temperature=100, # K, initial temperature 
                 optimizer=SDLBFGS,
                 fmax=0.05,
                 dr=0.4,
                 logfile='-', 
                 trajectory=None,
                 optimizer_logfile='-',
                 local_minima_trajectory='local_minima.con',
                 adjust_cm=True,
                 mss=0.1,
                 minenergy=None,
                 distribution='uniform',
                 adjust_step_size=None,
                 adjust_every =10,
                 target_ratio = 0.5,
                 adjust_fraction = 0.05,
                 significant_structure = True,  # displace from minimum at each move
# JD: removing this flag;    significant_structure2 = False, # displace from global minimum found so far at each move
                 pushapart = 0.4,
                 jumpmax=10,
                 jmp = 7, # number of jump steps taken in BHOJ
    		 molecular_dynamics = False, # trial move using md when true
		 dimer_a = 0.001,
		 dimer_d = 0.01,
		 dimer_steps = 20,
		 timestep = 0.1, # moleculare dynamics time step
		 mdmin = 2, # number of minima to pass in md before stopping
		 history_weight = 0.0, # the weight factor of history >= 0 
                 adjust_temp = True, # dynamically adjust the temperature in BH acceptance
		 minimaHopping_acceptance = False, # use MH acceptance criteria instead of BH
                 minimaHopping_history = True, # use history in MH acceptance criteria
		 beta1 = 1.04, # temperature adjustment parameter
		 beta2 = 1.04, # temperature adjustment parameter
		 beta3 = 1.0 / 1.04, # temperature adjustment parameter
	         Ediff0 = 0.5,  # eV, initial energy acceptance threshold
        	 alpha1 = 0.98,  # energy threshold adjustment parameter
        	 alpha2 = 1. / 0.98,  # energy threshold adjustment parameter
		 minima_threshold = 2,  # round energies to how many decimal place
                 use_geometry = False,
                 eps_r = 0.1,
                 confile_directory = "confiles" # cannot be a directory that exists
                 
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
#        self.significant_structure2 = significant_structure2
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
        self.use_geo = use_geometry
        self.eps_r = eps_r
        self.con_dir = confile_directory

        # when a MD sim. has passed a local minimum:
        self.passedminimum = PassedMinimum()

        self.mss = mss
	# dictionary for found local minima
        # keys will be the potential energy rounded to self.minima_threshold digits left of the decimal
        # values will be number of times the potential energy has been visited 
	self.minima = {}
        # dictionary for geometry comparison
        # keys = "#.con", values = # visits
        self.cons = {}
        self.initialize()

    def initialize(self):
        self.positions = 0.0 * self.atoms.get_positions()
        self.Emin = self.get_energy(self.atoms.get_positions()) or 1.e32 
        self.rmin = self.atoms.get_positions()
        self.positions = self.atoms.get_positions()
        self.call_observers()
        self.log(-1, self.Emin, self.Emin,self.dr)
        if self.use_geo:
            os.mkdir(self.con_dir)
                
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
        approxEn = round(En,self.minima_threshold)
        approxEo = round(Eo,self.minima_threshold)
        if approxEn in self.minima:
            self.minima[approxEn] += 1
            return False, approxEn, approxEo
        else:
            self.minima[approxEn] = 1
            return True, approxEn, approxEo

    def pbc(self, r, box, ibox = None):
        """
        Applies periodic boundary conditions.
        Parameters:
        r:      the vector the boundary conditions are applied to
        box:    the box that defines the boundary conditions
        ibox:   the inverse of the box. This will be calcluated if not provided.
        """
        if ibox is None:
            ibox = np.linalg.inv(box)
        vdir = np.dot(r, ibox)
        vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
        return np.dot(vdir, box)

    def per_atom_norm(self, v, box, ibox = None):
        '''
        Returns a length N numpy array containing per atom distance
        v:	an Nx3 numpy array
        box:    box matrix that defines the boundary conditions
        ibox:   the inverse of the box. will be calculated if not provided
        '''
        diff = self.pbc(v, box, ibox)
        return np.sqrt(np.sum(diff**2.0, axis=1))

    def rot_match(self, a, b, eps_r):
        """Taken from eon software package. Determines if
        atoms a is the same geometry as atoms b."""
        if not (a.free.all() and b.free.all()):
            print "Comparing structures with frozen atoms with rotational matching; check_rotation may be set incorrectly"
        acm = sum(a.r)/len(a)
        bcm = sum(b.r)/len(b)

        ta = a.copy()
        tb = b.copy()
        ta.r -= acm
        tb.r -= bcm

        #Horn, J. Opt. Soc. Am. A, 1987
        m = np.dot(tb.r.transpose(), ta.r)
        sxx = m[0][0]
        sxy = m[0][1]
        sxz = m[0][2]
        syx = m[1][0]
        syy = m[1][1]
        syz = m[1][2]
        szx = m[2][0]
        szy = m[2][1]
        szz = m[2][2]

        n = np.zeros((4,4))
        n[0][1] = syz-szy
        n[0][2] = szx-sxz
        n[0][3] = sxy-syx

        n[1][2] = sxy+syx
        n[1][3] = szx+sxz

        n[2][3] = syz + szy

        n += n.transpose()

        n[0][0] = sxx + syy + szz
        n[1][1] = sxx-syy-szz
        n[2][2] = -sxx + syy -szz
        n[3][3] = -sxx -syy + szz

        w,v = np.linalg.eig(n)
        maxw = 0
        maxv = 0
        for i in range(len(w)):
            if w[i] > maxw:
                maxw = w[i]
                maxv = v[:,i]

        R = np.zeros((3,3))

        aa = maxv[0]**2
        bb = maxv[1]**2
        cc = maxv[2]**2
        dd = maxv[3]**2
        ab = maxv[0]*maxv[1]
        ac = maxv[0]*maxv[2]
        ad = maxv[0]*maxv[3]
        bc = maxv[1]*maxv[2]
        bd = maxv[1]*maxv[3]
        cd = maxv[2]*maxv[3]

        R[0][0] = aa + bb - cc - dd
        R[0][1] = 2*(bc-ad)
        R[0][2] = 2*(bd+ac)
        R[1][0] = 2*(bc+ad)
        R[1][1] = aa - bb + cc - dd
        R[1][2] = 2*(cd-ab)
        R[2][0] = 2*(bd-ac)
        R[2][1] = 2*(cd+ab)
        R[2][2] = aa - bb - cc + dd
        tb.r = np.dot(tb.r, R.transpose())

        dist = max(self.per_atom_norm(ta.r - tb.r, ta.box))
        return dist < eps_r

    def find_match(self):
        """ determines if atoms is the same geometry as any previously visited minima."""
        new_con = str(len(self.cons.keys())) + ".con"
        tsase.io.write_con(new_con, self.atoms, w='w')
        a = fileio.loadcon(new_con)
        for k in self.cons.keys():
            b = fileio.loadcon(k)
            #same, r = geometry.identical(a, b, self.eps_r)
            #same = geometry.match(a, b, self.eps_r, 1.5, True, True, False)
            same = self.rot_match(a, b, self.eps_r)
            #print "r", r
            if same:
                os.remove(new_con)
                #print "match", k
                return k
        # atoms is unvisited local minima 
        os.remove(new_con)
        return None

    def update_con(self, confile):
        """ update dictionary of .con files
        add a new .con if atoms is a new local min
        or increment the # of visits for the correct .con """
        if confile is not None:
            # debugging
            localMin = tsase.io.read_con(confile)
            lj = tsase.calculators.lj(cutoff=35.0)
            localMin.center(100.0)
            localMin.set_calculator(lj)
            #print "repeat min:", self.atoms.get_potential_energy()
            #print "equivalent to", localMin.get_potential_energy()
            self.cons[confile] += 1
        else:
            new_con = self.con_dir + "/" + str(len(self.cons.keys())) + ".con"
            tsase.io.write_con(new_con, self.atoms, w='w')
            self.cons[new_con] = 1

    def _maxwellboltzmanndistribution(self,masses,N,temp,communicator=world):
        # For parallel GPAW simulations, the random velocities should be
        # distributed.  Uses gpaw world communicator as default, but allow
        # option of specifying other communicator (for ensemble runs)
        xi = np.random.standard_normal((len(masses), 3))
        if N is not None:
            xi = N
            #momenta = 3 * len(masses) * xi * np.sqrt(2 * masses * temp)[:, np.newaxis]
            #momenta = xi * np.sqrt((3/2) * len(masses) * masses * (self.temperature * kB))[:, np.newaxis]
        #else:
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
        
        #return self.energy
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
       # JD: Removed to be consistent with the corrected basin.py
       #     if self.significant_structure == True:
       #         rn = self.local_min_pos + disp
       #     elif self.significant_structure2 == True:
       #         ro = self.get_minimum()
       #         rn = ro + disp
       #     else:
            rn = ro + disp
            rn = self.push_apart(rn)
            self.atoms.set_positions(rn)
            if self.cm is not None:
                cm = self.atoms.get_center_of_mass()
                self.atoms.translate(self.cm - cm)
        else :
            # Move atoms using Dimer and an MD step
            dimer = ModifiedDimer()
            N = dimer(self.atoms, self.dimer_a, self.dimer_d, self.dimer_steps)
            self._molecular_dynamics(step, N)
        # JD: Add displacement to ro at each step regardless of significant structure flag
        rn = ro + disp
        rn = self.atoms.get_positions()
        world.broadcast(rn, 0)
        self.atoms.set_positions(rn)
        return self.atoms.get_positions()

    def acceptance_MH(self, steps, ro, Eo, maxtemp):
        """Adjusts parameters and positions based  on minima hopping acceptance criteria."""
        rejectnum = 0
        lastcon = "none"
        for step in range(steps):
            positionsOld = self.atoms.get_positions()
            En = None
            rn = self.move(step,ro)
            En = self.get_energy(rn)
            self.steps += 1
            while En is None:
                rn = self.move(step,ro)
                self.atoms.set_positions(rn)
                En = self.get_energy(rn)
            if En < self.Emin:
                self.Emin = En
                self.rmin = self.atoms.get_positions()
                self.call_observers()
            self.log(step, En, self.Emin,self.dr)
            match = None
            countEo = 0
            countEn = 0
            if self.use_geo:
                match = self.find_match()
            else:
                match, countEn, countEo = self.find_ene_match(En, Eo)
            if match is not None:
                if (self.use_geo and match == lastcon) or match == 1:
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
            if self.jumpmax and (rejectnum > self.jumpmax):
                #JMP???
                for i in range(0,self.jmp):
                    rn = self.move(step,rn)
                accept = True
            if accept:
                rejectnum = 0
                ro = rn
                if self.use_geo:
                    new_con = self.update_con(match)
                    lastcon = new_con
                else:
                    self.update_minima(En, Eo)
                Eo = En
                if self.lm_trajectory is not None:
                    tsase.io.write_con(self.lm_trajectory,self.atoms,w='a')
            else:
                self.atoms.set_positions(positionsOld)
            if self.minenergy is not None:
                if Eo < self.minenergy:
                    #print "geo: ", self.cons.values()
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
            lastcon = "none"
            while En is None:
                rn = self.move(step,ro)
                En = self.get_energy(rn)
            if En < self.Emin:
                self.Emin = En
                self.rmin = self.atoms.get_positions()
                self.call_observers()
            self.log(step, En, self.Emin,self.dr)
            match = None
            countEo = 0
            countEn = 0
            if self.use_geo:
                match = self.find_match()
            else:
                match, countEn, countEo = self.find_ene_match(En, Eo)
            if self.adjust_temp and match is not None:
                if (self.use_geo and match == lastcon) or match == 1:
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
                if match is not None:
                    if self.use_geo:
                        hn = self.cons[match] / step
                    else:
                        hn = countEn / step
                if self.use_geo and lastcon != "none":
                    ho = self.cons[lastcon] / step
                else:
                    ho = countEo / steps
                kT = self.temperature * kB
                # occasionally overflowing the floating point :(
                #accept = np.exp(((Eo - En) + (self.w * (ho - hn))) / kT) > np.random.uniform()
                val = ((Eo - En) + (self.w * (ho - hn))) / kT
                if val < 1.0:
                    accept = np.exp(val) > np.random.uniform()
                else: # accept the new position
                    accept = True;
            if self.jumpmax and rejectnum > self.jumpmax:
                #JMP???
                for i in range(0,self.jmp):
                    rn = self.move(step,rn)
                accept = True
            if accept:
                acceptnum += 1.
                recentaccept += 1.
                rejectnum = 0
                ## JD: made edit to set ro to local minima if we accept the configuration.  eliminating significant_structure 2 flag
                if self.significant_structure == True:
                    ro = self.local_min_pos.copy()
                else:
                    ro = rn.copy()
                if self.use_geo:
                    new_con = self.update_con(match)
                    lastcon = new_con
                else:
                    self.update_minima(En, Eo)
                Eo = En
                if self.lm_trajectory is not None:
                    tsase.io.write_con(self.lm_trajectory,self.atoms,w='a')
            else:
                rejectnum += 1
                self.atoms.set_positions(positionsOld)
            if self.minenergy != None:
                if Eo < self.minenergy:
                    #print "geo: ", self.cons.values()
                    break
            if self.adjust_step == True:
                if step % self.adjust_every == 0:
                    #ratio = float(acceptnum)/float(self.steps)
                    ratio = float(acceptnum)/float(steps)
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
        print "geo: ", self.cons.values()
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

