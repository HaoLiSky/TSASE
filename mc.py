#!/usr/bin/env python
import numpy
from ase import Atoms
from tsase.io import read_con
from ase import units

class MonteCarlo:
    def __init__(self, atoms, temperature):
        self.atoms = atoms
        self.callback_functions = []
        self.beta = 1.0/(units.kB * temperature)

        self.steps = 0
        self.accepts = 0
        self.rejects = 0

    def get_acceptance_ratio(self):
        return float(self.accepts)/float(self.steps)

    def attach(self, callback, interval=1):
        self.callback_functions.append((callback, interval))

    def run(self, steps, displacement_size=0.1):
        self.steps = 0
        self.accepts = 0
        self.rejects = 0
        self.e = self.atoms.get_potential_energy()
        for i in xrange(steps):
            r = self.atoms.get_positions()
            dr = numpy.random.normal(0.0, displacement_size, r.shape)
            r_trial = r + dr
            self.atoms.set_positions(r_trial)
            e_trial = self.atoms.get_potential_energy()
            self.atoms.set_positions(r)

            delta_e = e_trial - self.e

            accept = False
            if delta_e < 0.0:
                accept = True
            else:
                U = numpy.random.random()
                if U < numpy.exp(-self.beta*delta_e):
                    accept = True

            if accept:
                self.accepts += 1
                self.atoms.set_positions(r_trial) 
                self.e = e_trial
            else:
                self.rejects += 1

            self.steps += 1

            for callback in self.callback_functions:
                function, interval = callback
                if self.steps % interval == 0:
                    function(self)
