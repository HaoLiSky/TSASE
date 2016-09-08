#!/usr/bin/env python

import numpy

from expectra.expectra import Expectra
from expectra.exafs import exafs_first_shell, exafs_multiple_scattering
from expectra.io import read_xdatcar, read_con, read_chi
from expectra.feff import load_chi_dat
from ase import atoms

def main():
    filename = 'CONTCAR'
    trajectory = read_con(filename)
    calc = expectra(atoms=trajectory, kmax = 10.0)
    trajectory.set_calculator(calc)
    chi_devia = trajectory.get_chi_deviation()
    print(chi_devia)
if __name__ == '__main__':
    main()
    
