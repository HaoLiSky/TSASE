#!/usr/bin/env python
import sys
sys.path.append("..")
import os

from ase.test import NotAvailable
from ase.io import read
from ase.optimize.lbfgs import LBFGS
from ase.calculators.emt import EMT
from expectra.aselite import Atoms
from expectra.cal_exafs import Expectra
from expectra.io import read_xdatcar, read_con, read_chi
#from ase import *
from expectra.basin import BasinHopping

def main():
    p1 = read('geometry.xyz',index=0,format='xyz')
    p1.set_cell([[20,0,0],[0,20,0],[0,0,20]],scale_atoms=False,fix=None)
    p1.set_pbc((True, True, True))
#    filename='CONTCAR'
#    p1 = read_xdatcar(filename)
#    print(p1[0].get_positions())
#    print "p1 type in test.py"
#    print(type(p1))
    exafs_calc = Expectra(kmax = 10.0)
    bh = BasinHopping(atoms=p1,
                      opt_calculator = EMT(),
                      exafs_calculator=exafs_calc,
#                      temperature=300 * kB,
                      dr=0.5,
                      logfile='pot_log',
                      optimizer=LBFGS,
                      fmax=0.05)
    bh.run(2)
#    p1[0].set_calculator(exafs_calc)
#    chi_devia = p1[0].get_potential_energy()
#    print(chi_devia)
if __name__ == '__main__':
    main()
    
