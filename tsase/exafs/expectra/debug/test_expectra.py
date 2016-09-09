#!/usr/bin/env python

from expectra.cal_exafs import Expectra
from expectra.io import read_xdatcar, read_con, read_chi
#from ase import atoms
from ase.io.vasp import read_vasp

def main():
    filename='CONTCAR'
    p1 = read_xdatcar(filename)
#    print(p1[0].get_positions())
#    print(type(p1[0]))
    calc = Expectra(atoms=p1[0], kmax = 10.0)
    p1[0].set_calculator(calc)
    chi_devia = p1[0].get_potential_energy()
    print(chi_devia)
if __name__ == '__main__':
    main()
    
