#!/usr/bin/env python
'''
New calculator by xph at 05/09/2016
applying external force to any existing calculator
'''

import numpy as np

class applyF():
    def __init__(self, calc, pullatoms=None, pullF=None):
        """
        calc:      the original calculator
        pullatoms: a list of atom identfication numbers to be pulled.
                   For example, [0, -1] is to pull the first and last atoms.
        pullF:     pulling forces, a numpy array with the shape of (len(pullatoms), 3). 
                   There should be: len(pullF) = len(pullatoms);
                   Each component of pullF is a 3D vector.
        """
        if pullatoms != None:
            self.pullatoms = pullatoms
            self.pullF     = pullF
            if pullF.shape != (len(pullatoms), 3):
                raise ValueError('pullF must be in the shape of (len(pullatoms),3)')
            elif np.sum(pullF) > 0.0001:
                print np.sum(pullF)
                raise ValueError('Applied force is not balanced: summation not zero')
        else:
            raise ValueError('Which atom to pull?')

        self.calc = calc

    # Pipe all the stuff from calc that is not overwritten.
    # Pipe all requests for get_original_* to self.calc.
    def __getattr__(self, attr):
        """Return any value of the Atoms object"""
        return getattr(self.calc, attr)

    def get_forces(self, atoms):
        f    = self.calc.get_forces(atoms)
        for i in range(len(self.pullatoms)):
            atomnum = self.pullatoms[i]
            f[atomnum] += self.pullF[i]
        return f
            
    def get_potential_energy(self, atoms):
        e    = self.calc.get_potential_energy(atoms)
        r    = atoms.positions[np.array(self.pullatoms)]
        e   -= np.vdot(self.pullF, r)
        #for i in range(len(self.pullatoms)):
        #    atomnum = self.pullatoms[i]
        #    ri      = atoms.positions[atomnum]
        return e

    
