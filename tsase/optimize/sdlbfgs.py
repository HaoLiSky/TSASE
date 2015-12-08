# -*- coding: utf-8 -*-
import sys
import numpy as np
from tsase.optimize.optimize import Optimizer

class SDLBFGS(Optimizer):
    """Limited memory BFGS optimizer.
    
    A limited memory version of the bfgs algorithm. Unlike the bfgs algorithm
    used in bfgs.py, the inverse of Hessian matrix is updated.  The inverse
    Hessian is represented only as a diagonal matrix to save memory

    This version of LBFGS is based off of ASE implementation with a few improvements 
    """
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=0.2, memory=100, damping = 1.0):
        """
        Parameters:

        restart: string
            Pickle file used to store vectors for updating the inverse of Hessian
            matrix. If set, file with such a name will be searched and information
            stored will be used, if the file exists.

        logfile: string
            Where should output go. None for no output, '-' for stdout.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        maxstep: float
            This parameter has been changed from ASEs implementation to a max total stepsize
            Default is 0.2 Angstrom.

        memory: int
            Number of steps to be stored. Default value is 100. Three numpy
            arrays of this length containing floats are stored.

        damping: float
            The calculated step is multiplied with this number before added to
            the positions. 
            
        """
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

        if maxstep is not None:
            if maxstep > 1.0:
                raise ValueError('You are using a much too large value for ' +
                                 'the maximum step size: %.1f Angstrom' % maxstep)
            self.maxstep = maxstep
        else:
            self.maxstep = 0.2

        self.memory = memory
        self.damping = damping
        self.p = None
        self.function_calls = 0
        self.force_calls = 0
        self.prev_force = 0
        self.prev_positions = 0

    def initialize(self):
        """Initalize everything so no checks have to be done in step"""
        self.iteration = 0
        self.s = []
        self.y = []
        self.rho = [] # Store also rho, to avoid calculationg the dot product
                      # again and again

        self.r0 = None
        self.f0 = None
        self.e0 = None
        self.task = 'START'
        self.load_restart = False

    def read(self):
        """Load saved arrays to reconstruct the Hessian"""
        self.iteration, self.s, self.y, self.rho, \
        self.r0, self.f0, self.e0, self.task = self.load()
        self.load_restart = True

    def step(self, f):
        """Take a single step
        Use the given forces, update the history and calculate the next step --
        then take it"""
        r = self.atoms.get_positions()
        p0 = self.p
        if self.iteration == 0:
            self.force_calls += 2
            ### initial Hessian by steepest descent step
            C = self.get_curvature_via_FD()
            H0 = 1/C
        else:
        ### change H0 to be approximate curvature at each step 
            self.force_calls += 1
            deltaF = self.prev_force - f 
            deltaR = r - self.prev_positions
            C = np.vdot(deltaF,deltaF)/np.vdot(deltaR,deltaF)
        ### if curvature is negative restart build up of hessian
        if C < 0:
            self.rho = []
            self.y = []
            self.s = []
        else:
            H0 = 1/C
            self.prev_force = f
            self.prev_positions = r
            self.update(r, f, self.r0, self.f0)
        s = self.s
        y = self.y
        rho = self.rho
        loopmax = len(y)
        #### if we reset hessian because of negative curvature take maxstep in the direction of the force	
        if C < 0:
            g = f*1000
            dr = self.determine_step(g) * self.damping
            self.atoms.set_positions(r+dr)
        ## otherwise take lbfgs step
        else:
            a = np.empty((loopmax,), dtype=np.float64)
            ### The algorithm itself:
            q = - f.reshape(-1) 
            for i in range(loopmax - 1, -1, -1):
                a[i] = rho[i] * np.dot(s[i], q)
                q -= a[i] * y[i]
            z = H0 * q
        
            for i in range(loopmax):
                b = rho[i] * np.dot(y[i], z)
                z += s[i] * (a[i] - b)

            #### check if angle of move is greater than 90 reset build up of Hessian #######
            stepdir = -z/np.sqrt(np.vdot(z,z))
            fdir = f/np.sqrt(np.vdot(f,f))
            dot = np.vdot(fdir,stepdir)
            if dot > 1:
                dot = 1
            if dot < -1:
                dot = -1
            angle = np.arccos(dot)*360/(2*np.pi)
            angle = 0.0 ######### added for test!!!!!!
            if angle > 90:  ### if greater than 90 take SD step
                self.rho = []
                self.y = []
                self.s = []
                g = f*1000 
                dr = self.determine_step(g) * self.damping
                self.atoms.set_positions(r+dr)
            else:  ##### otherwise take lbfgs step
                self.p = - z.reshape((-1, 3))
                g = -f
                dr = self.determine_step(self.p) * self.damping
                self.atoms.set_positions(r+dr)
 
        self.iteration += 1
        self.r0 = r
        self.f0 = -g
        self.dump((self.iteration, self.s, self.y, 
                   self.rho, self.r0, self.f0, self.e0, self.task))

    def get_curvature_via_FD(self):
        dR = 1.0e-5
        pos = self.atoms.get_positions()
        initforce = self.atoms.get_forces()
        Fnorm = initforce/np.sqrt(np.vdot(initforce,initforce))
        new_pos = pos + dR*Fnorm
        self.atoms.set_positions(new_pos)
        newforce = self.atoms.get_forces()
        C = -np.vdot(newforce - initforce,Fnorm)/dR
        self.atoms.set_positions(pos)
        return C

    def determine_step(self, dr):
        """Determine step to take according to maxstep

        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """
        steplengths = (dr**2).sum(1)**0.5
        longest_step = np.max(steplengths)
        if longest_step >= self.maxstep:
            dr *= self.maxstep / longest_step

        return dr


#    def determine_step(self, dr):
#        """Determine step to take according to maxstep
#        
#        Normalize all steps as the largest step. This way
#        we still move along the eigendirection.
#        """
#       length = np.sqrt(np.vdot(dr,dr))
#       if length > self.maxstep:
#           dr *= self.maxstep/length
#        
#        return dr

    def update(self, r, f, r0, f0):
        """Update everything that is kept in memory

        This function is mostly here to allow for replay_trajectory.
        """
        if self.iteration > 0:
            s0 = r.reshape(-1) - r0.reshape(-1)
            self.s.append(s0)

            # We use the gradient which is minus the force!
            y0 = f0.reshape(-1) - f.reshape(-1)
            self.y.append(y0)
            
            rho0 = 1.0 / np.dot(y0, s0)
            self.rho.append(rho0)

        if len(self.s) > self.memory:
            self.s.pop(0)
            self.y.pop(0)
            self.rho.pop(0)


    def replay_trajectory(self, traj):
        """Initialize history from old trajectory."""
        if isinstance(traj, str):
            from ase.io.trajectory import PickleTrajectory
            traj = PickleTrajectory(traj, 'r')
        r0 = None
        f0 = None
        # The last element is not added, as we get that for free when taking
        # the first qn-step after the replay
        for i in range(0, len(traj) - 1):
            r = traj[i].get_positions()
            f = traj[i].get_forces()
            self.update(r, f, r0, f0)
            r0 = r.copy()
            f0 = f.copy()
            self.iteration += 1
        self.r0 = r0
        self.f0 = f0

    def func(self, x):
        """Objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.function_calls += 1
        return self.atoms.get_potential_energy()

    def fprime(self, x):
        """Gradient of the objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.force_calls += 1
        # Remember that forces are minus the gradient!
        return - self.atoms.get_forces().reshape(-1)

