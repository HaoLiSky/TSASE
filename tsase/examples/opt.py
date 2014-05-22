#!/usr/bin/env python
'''

Cell optimization example with a shear of 3.0 GPa applied on zx direction
'''

from ase.optimize.fire import FIRE
from ase import *
from ase.io import read,write
import os
import sys
import numpy as np
from tsase.mushybox import mushybox
from tsase.calculators.vasp_ext import Vasp

p1 = read('POSCAR',format='vasp')

calc = Vasp(prec = 'Normal', 
            ediff = 1e-5,
            kpts = (4,3,4),
            gamma= True,
            voskown = 1,
            lcharg = False,
            isym = 0,
            ispin = 2,
            npar = 2,
            nsim = 4,
            algo = 'Fast',
            lreal= 'Auto',
            lplane = True,
            encut= 400,
            ismear = 0,
            sigma  = 0.08,
            ldau      = True,
            ldautype  = 2,
            ldauprint = 2,
            ldau_luj = {'Ca':{'L':-1,'U':0,'J':0},'Ir':{'L':2, 'U':3.75, 'J':1.0},'O':{'L':-1,'U':0,'J':0}},
            lmaxmix   = 4
              )
p1.set_calculator(calc)

pstress = p1.get_cell()*0.0
pstress[2][0] = 3.0

p1box = mushybox(p1,pstress)
print len(p1box)
print p1box.jacobian
print p1box.get_potential_energy()

dyn = FIRE(p1box)
#dyn = MDMin(p1box)
dyn.run(fmax=0.01)

io.write("CONTCAR",p1,format='vasp')

