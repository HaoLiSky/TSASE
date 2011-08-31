#!/usr/bin/env python

'''
ssneb optimization example
output: 1. fe.out
                  contais force, energy and distance between images
                  the final results can be plotted by nebspline.pl in vtstscripts
        2. i.CON: i is the number of image
                  geomety of each image, in vasp POSCAR format
'''

from tsase import neb
from tsase.calculators.lammps_ext import LAMMPS
from ase.io import read
import os

# read geometry and set charges
p1 = read('CdSe_hex',format='vasp')
p2 = read('CdSe_sq',format='vasp')
tags = [a.symbol == 'Se' for a in p1]
charges = [(-1)**i*1.18 for i in tags]
p1.set_charges(charges)
p2.set_charges(charges)

# set calculator to LAMMPS in tsase. ASE's LAMMPS calculator doesn't write charge information into data file.
pair_coeff = [ '1 1 0.00145 1.98', '2 2 0.00128 5.24' ]
parameters = { 'pair_style':'lj/cut/coul/long 10.0 10.0', 'pair_coeff':pair_coeff, 'kspace_style':'ewald/n 1.0e-8', 'atom_style':'charge','mass':['1 1','2 1'], 'pair_modify':'table 12 mix arithmetic'}
calc = LAMMPS(parameters=parameters)
# or set calculator to VASP
#calc = Vasp(prec = 'Accurate',
#            ediff = 1e-6,
#            kpts = (6,6,6),
#            voskown = 1,
#            lcharg = False,
#            algo = 'Fast',
#            lreal= False,
#            lplane=True
#              )

p1.set_calculator(calc)
p2.set_calculator(calc)

# initialize the band 
nim = 7  # number of images, including end points
band = neb.ssneb(p1, p2, numImages = nim, method = 'ci', weight = 1)
# to restart uncomment the following lines, which read the previous optimized images into the band
#for i in range(1,nim-1):
#    filename = str(i)+'.CON'
#    b = read(filename,format='vasp')
#    band.path[i].set_positions(b.get_positions())
#    band.path[i].set_cell(b.get_cell())


#opt = neb.qm_ssneb(band, maxmove = 0.20, dt = 0.05)
opt = neb.fire_ssneb(band, maxmove =0.2, dtmax = 0.1, dt=0.1)
opt.minimize(forceConverged=0.001, maxIterations = 1000)



