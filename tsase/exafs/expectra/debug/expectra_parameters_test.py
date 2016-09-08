#!/usr/bin/env python

#from __future__ import print_function
import mpi4py.MPI

import argparse
import sys
import numpy

from expectra.exafs import exafs_first_shell, exafs_multiple_scattering
from expectra.io import read_xdatcar, read_con, read_chi
from expectra.feff import load_chi_dat
from ase.calculators.calculator import Calculator, all_changes, Parameters


"""
Calculator is the superclass or base class. expectra is the subclass or derived
class
"""

class expectra(Calculator):

    implemented_properties = ['chi_deviation']

    default_parameters = dict(
        neighbor_cutoff = 6.0,
        s02 = 0.89,
        energy_shift = 3.4,
        edge = 'L3',
        absorber = 'Au',
        skip = 0,
        exp_chi_fie = 'chi_exp.dat')

    def __int__(self, label='EXAFS', 
                atoms=None, kmin=0.00, kmax=10.00, chi_deviation=100, **kwargs):
                
        """
        The expectra constructor:
        kmin, kmax ...........define the k window you are interested on. It is
        suggested to set them as the values appeared in experimental data
        atoms.................coordinates or trajectories
        chi_deviaiton.........used to store the calcuation deviation between
        experimental and calculated EXAFS spectra
        """
        self.label = label
        self.atoms = atoms
        self.kmin = kmin
        self.kmax = kmax
        self.chi_deviation = chi_deviation
        self.parameters = None
        self.results = None
    
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms,
                            **kwargs)
#        print(dict.get(neighbor_cutoff))
#        print(dict.get(s02))
#        print(dict.get(energy_shift))


    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2,
        ...)."""
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
           self.reset()

def main():
    test = expectra(kmin = 2.00,
                    kmax = 10.00,
                    chi_deviation = 10,
                    neighbor_cutoff = 2.0,
                    s02 = 0.5)
    f = open('out.dat', 'w')
    f.write("%6.3f\n" % (test.parameters.neighbor_cutoff))
    f.write("%6.3f\n" % (test.parameters.chi_deviation))
    f.write("%6.3f\n" % (test.parameters.s02))
    f.close()

if __name__ == '__main__':
    main()
