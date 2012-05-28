#!/usr/bin/env python
'''

Cell optimization class
'''

from ase.optimize.fire import FIRE
from ase import *
from ase.io import read,write
from ase.calculators.lammps import LAMMPS
import os
import sys
import numpy as np

class mushybox:
    def __init__(self, atomsx, express=np.zeros((3,3))):
        """box relaxation
        atomsx: an Atoms object
        express: external pressure, a 3*3 lower triangular matrix in the unit of GPa
                 define positive values as compressing
        """
        self.atomsx = atomsx 
        self.express= express * units.GPa
        if express[0][1]**2+express[0][2]**2+express[1][2]**2 > 1e-5:
           print "warning: xy, xz, yz components of the external pressure will be set to zero"

        cell       = atomsx.get_cell()
        vol        = atomsx.get_volume()
        self.natom = atomsx.get_number_of_atoms()
        avglen     = (vol/self.natom)**(1.0/3.0)
        self.jacobian = avglen * self.natom**0.5

    def get_positions(self):
        r    = self.atomsx.get_positions()
        Rc   = np.vstack((r, self.atomsx.get_cell()*0.0))
        return Rc

    def set_positions(self,r):
        ratom  = r[:-3]
        self.atomsx.set_positions(ratom)
        rcell  = self.atomsx.get_cell()
        rcell += np.dot(rcell,r[-3:])/self.jacobian
        self.atomsx.set_cell(rcell,scale_atoms=True)
        print "set_positions"
    
    def __len__(self):
        return self.natom+3

    def get_forces(self):
        #f    = self.atomsx.get_forces(apply_constraint)
        f    = self.atomsx.get_forces()
        stt  = self.atomsx.get_stress()
        vol  = self.atomsx.get_volume()*(-1)
        st   = np.zeros((3,3))
        #following the order of get_stress in vasp.py
        st[0][0] = stt[0] * vol  
        st[1][1] = stt[1] * vol
        st[2][2] = stt[2] * vol
        st[2][1] = stt[3] * vol
        st[2][0] = stt[4] * vol
        st[1][0] = stt[5] * vol
        st  -= self.express * (-1)*vol
        print "get_forces"
        print st
        Fc   = np.vstack((f, st/self.jacobian))
        return Fc
    
    def get_potential_energy(self):
        return self.atomsx.get_potential_energy()


