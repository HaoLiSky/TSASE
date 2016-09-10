#!/usr/bin/env python
import sys
sys.path.append("..")

from expectra.cal_exafs import Expectra
from expectra.io import read_xdatcar, read_con, read_chi
#from ase import atoms
from ase import *
from expectra.basin import BasinHopping

def main():
    filename='CONTCAR'
    p1 = read_xdatcar(filename)
#    print(p1[0].get_positions())
#    print(type(p1[0]))
    exafs_calc = Expectra(kmax = 10.0)
    bh = BasinHopping(atoms=p1[0],
                      exafs_calculator=exafs_calc,
#                      temperature=300 * kB,
                      dr=0.5,
#                      optimizer=LBFGS,
                      fmax=0.1)
#    p1[0].set_calculator(exafs_calc)
#    chi_devia = p1[0].get_potential_energy()
#    print(chi_devia)
if __name__ == '__main__':
    main()
    
