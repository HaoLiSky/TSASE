#!/usr/bin/env python

from expectra.cal_exafs import Expectra
from expectra.io import read_xdatcar, read_con, read_chi
from ase import atoms
from ase.io.vasp import read_vasp

def main():
    trajectory = read_vasp(filename='CONTCAR')
    calc = Expectra(atoms=trajectory, kmax = 10.0)
    trajectory.set_calculator(calc)
    chi_devia = trajectory.get_potential_energy()
    print(chi_devia)
if __name__ == '__main__':
    main()
    
