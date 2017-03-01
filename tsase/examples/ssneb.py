#!/usr/bin/env python

'''
ssneb optimization example
output: 1. fe.out
                  contais force, energy and distance between images
                  the final results can be plotted by nebspline.pl in vtstscripts
        2. i.CON: i is the number of image
                  geomety of each image, in vasp POSCAR format

If parallel=True is set, the ssneb calculation is parallelized over images. 
The command to execute this script needs to be:
    mpirun -np 5 python ssneb.py
where 5 is the number of intermidiate images (nim-2)
'''

from tsase import neb
from tsase.calculators.lammps_ext import LAMMPS
#from tsase.calculators.vasp_ext import Vasp
from ase.io import read
import os, numpy

'''
#-------------- uncomment if parallelize over image ------------------
# parallelize over images with mpi4py
from mpi4py import MPI
'''

# read geometry and set charges
p1 = read('CdSe_hex',format='vasp')
p2 = read('CdSe_sq',format='vasp')
tags = [a.symbol == 'Se' for a in p1]
#charges = [(-1)**i*1.18 for i in tags]
p1.set_tags(tags)
p2.set_tags(tags)
for p in p1: 
    p.charge = (-1)**p.tag*1.18
for p in p2: 
    p.charge = (-1)**p.tag*1.18

# set calculator to LAMMPS in tsase. ASE's LAMMPS calculator doesn't write charge information into data file.
pair_coeff = [ '1 1 0.00145 1.98', '2 2 0.00128 5.24' ]
parameters = { 'pair_style':'lj/cut/coul/long 10.0 10.0', 'pair_coeff':pair_coeff, 'kspace_style':'ewald 1.0e-8', 'atom_style':'charge','mass':['1 1','2 1'], 'pair_modify':'table 12 mix arithmetic'}
calc = LAMMPS(parameters=parameters, tmp_dir='trash')

# or set calculator to VASP
#calc = Vasp(prec = 'Normal',
#            ediff = 1e-5,
#            kpts = (6,6,6),
#            voskown = 1,
#            lcharg = False,
#            algo = 'Fast',
#            lreal= False,
#            lplane=True
#              )

p1.set_calculator(calc)
p2.set_calculator(calc)

# external stress applied in the unit of GPa
stress=numpy.zeros((3,3))
#stress[2,2] = -1.0  #negative is tension

# initialize the band 
nim = 7  # number of images, including end points
# no climbing image first
band = neb.ssneb(p1, p2, numImages = nim, express = stress)
#band = neb.ssneb(p1, p2, numImages = nim, method = 'ci', express = stress)
'''
#-------------- substitute for the above line if parallelize over image ------------------
band = neb.ssneb(p1, p2, numImages = nim, method = 'ci', express = stress, parallel = True)
'''

# to restart, uncomment the following lines which read the previous optimized images into the band
#for i in range(1,nim-1):
#    filename = str(i)+'.CON'
#    b = read(filename,format='vasp')
#    band.path[i].set_positions(b.get_positions())
#    band.path[i].set_cell(b.get_cell())


#opt = neb.qm_ssneb(band, maxmove = 0.10, dt = 0.05)
opt = neb.fire_ssneb(band, maxmove =0.1, dtmax = 0.1, dt=0.1)
opt.minimize(forceConverged=0.01, maxIterations = 300)

