#!/usr/bin/env python

'''
SSDimer example
'''

from tsase.dimer import ssdimer
from tsase.neb.util import vunit
from tsase.calculators.lammps_ext import LAMMPS
from ase.io import read, write
import os
import sys
import numpy as np

p = read('CdSe_hex',format='vasp')

natom = len(p)
vol   = p.get_volume()
jacob = (vol/natom)**(1.0/3.0) * natom**0.5

pfin = read('CdSe_sq',format='vasp')
mode = np.zeros((len(p)+3,3))
mode[-3:] = pfin.get_cell()-p.get_cell()
icell= np.linalg.inv(p.get_cell())
mode[-3:] = np.dot(icell, mode[-3:]) * jacob

mode = vunit(mode)
cellt = p.get_cell()+np.dot(p.get_cell(), mode[-3:]/jacob)
p.set_cell(cellt, scale_atoms=True)


tags = [a.symbol == 'Se' for a in p]
charges = [(-1)**i*1.18 for i in tags]
p.set_charges(charges)

pair_coeff = [ '1 1 0.00145 1.98', '2 2 0.00128 5.24' ]
parameters = { 'pair_style':'lj/cut/coul/long 10.0 10.0', 'pair_coeff':pair_coeff, 'kspace_style':'ewald/n 1.0e-8', 'atom_style':'charge','mass':['1 1','2 1'], 'pair_modify':'table 12 mix arithmetic'}
calc = LAMMPS(parameters=parameters)
p.set_calculator(calc)

d = ssdimer.SSDimer_atoms(p, mode = mode, rotationMax = 4, phi_tol=30)
d.search(minForce = 0.0001, movie = "dimer2.movie", interval = 20 )
write("dimer1.con", d.R0, format='vasp')
print p.get_potential_energy()
