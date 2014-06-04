"""Andersen NVT dynamics class."""

import sys, ase
import numpy as np
from ase.md.md import MolecularDynamics
from ase.parallel import world

class nvtandersen(MolecularDynamics):

    def __init__(self,atoms,timestep,temperature,alpha = 0.8,tcol=0.35,fixcm = True,
        trajectory=None,logfile=None,loginterval=1, communicator=world,hyperplane=False,eig=None):

        MolecularDynamics.__init__(self,atoms,timestep,trajectory,logfile,loginterval)

        self.alpha = alpha
        self.temperature = temperature # note that this is in units of kB
        self.fixcm = fixcm
        self.communicator = communicator
        self.dt = timestep
        self.tcol = tcol
        self.hyperplane = hyperplane
        self.eig = eig


    def set_temperature(self,temperature):
        self.temperature = temperature

    def get_temperature(self):
        return self.temperature

    def set_timestep(self,timestep):
        self.dt = timestep

    def get_timestep(self):
        return self.dt

    def vrand(self,v):
        """
        Returns a random vector with the same magnitude and shape as v
        """
        vtemp = np.random.randn(v.size)
        return vtemp.reshape(v.shape)

    def apply_thermostat(self):
        alpha = self.alpha           # collision strength
        dt = self.dt
        atoms = self.atoms
        kT = self.temperature  
        Tcol = self.tcol            # avg time between collisions
        Pcol = 1.0 - np.exp(-dt/Tcol)  # Probability of collisions
        n=0
        # temp vector to reset velocity
        masses = self.atoms.get_masses()
        if not self.atoms.has('momenta'):
            from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
            MaxwellBoltzmannDistribution(self.atoms, self.temperature)
        tmpv = atoms.get_velocities() 
        for j in tmpv:
            if (np.random.random_sample() < Pcol): # check for collision
                new_v = np.sqrt(kT/masses[n]) * self.vrand(j)  # call tsse.util for Gaussian rand vector
                tmpv[n] = np.sqrt(1 - alpha ** 2) * j + alpha * new_v # mix old and rand velocities
            n += 1
        if self.hyperplane == True:
            tmpv -= np.vdot(tmpv,self.eig)*self.eig    
        self.atoms.set_velocities(tmpv)
        return

    def step(self, f):

        self.apply_thermostat()

        atoms = self.atoms
        p = self.atoms.get_momenta()
        p += 0.5 * self.dt * f

        if self.fixcm:
            psum = p.sum(axis=0) / float(len(p))
            p = p - psum

        self.atoms.set_positions(self.atoms.get_positions() + self.dt * p / self.atoms.get_masses()[:,np.newaxis])
        self.atoms.set_momenta(p)
        f = self.atoms.get_forces()
        atoms.set_momenta(self.atoms.get_momenta() + 0.5 * self.dt * f)

        return f
