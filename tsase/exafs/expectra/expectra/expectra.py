#!/usr/bin/env python
import mpi4py.MPI

import argparse
import sys

import numpy

from expectra.exafs import exafs_first_shell, exafs_multiple_scattering
from expectra.io import read_xdatcar, read_con

COMM_WORLD = mpi4py.MPI.COMM_WORLD

def mpiexcepthook(type, value, traceback):
    sys.__excepthook__(type, value, traceback)
    sys.stderr.write("exception occured on rank %i\n" % COMM_WORLD.rank)
    COMM_WORLD.Abort()
sys.excepthook = mpiexcepthook

class expectra:
    def __int__(self, kmin, kmax, trajectory, neighbor_cutoff = 6.0, s02 =
                0.89, energy_shift = 3.4, edge = 'L3', absorber = 'Au', skip =
                0, every = 1,multiple_scattering = True, exp_chi_file='chi_exp.dat')
    """
    The expectra constructor:
        kmin, kmax ...........define the k window you are interested on. It is
        suggested to set them as the value appeared in experimental data
        trajectory......coordinates or trajectories
    """
    self.kmin = kmin
    self.kmax = kmax
    self.trajectory = trajectory
    self.neighbor_cutoff = neighbor_cutoff
    self.s02 = s02
    self.energy_shift =energy_shift
    self.edge = edge
    self.absorber = absorber
    self.skip = skip
    self.every = every

#load experimental chi data
    try:
        k_exp, chi_exp = read_chi(exp_chi_file) 
    except:
        k_exp, chi_exp = load_chi_dat(exp_chi_file)

#    last_index = len(k_exp) - 1
#    kmin = k_exp[0]
#    kmax = k_exp[last_index]
    dk   = k_exp[1] - k_exp[0]

#    if args.ignore_elements:
#        args.ignore_elements = args.ignore_elements.split(',')

    #calculate experiment data
    trajectory = COMM_WORLD.bcast(self.trajectory)

    self.absorber = get_default_absorber(trajectory[0], self)

    k, chi = exafs_trajectory(self, trajectory)

    save_result(k, chi)

    #modify it to save a history of (chi, k)
    def save_result(k, chi):
        if COMM_WORLD.rank != 0: return
        print 'saving result to chi.dat'
        f = open('chi.dat', 'w')
        for i in xrange(len(k)):
            f.write("%6.3f %16.8e\n" % (k[i], chi[i]))
        f.close()
    
    def get_default_absorber(atoms, args):
        symbols = set(atoms.get_chemical_symbols())
        if args.absorber:
            if args.absorber not in symbols:
                print 'ERROR: --absorber %s is not in the system' % args.absorber
                sys.exit(2)
            else:
                return args.absorber
        if args.ignore_elements:
            symbols -= set(args.ignore_elements)
        if len(symbols) == 1:
            return list(symbols)[0]
        else:
            print 'ERROR: must specify --absorber if more than one chemical specie'
            sys.exit(2)
 
    def exafs_trajectory(args, trajectory):
        if args.multiple_scattering:
            k, chi = exafs_multiple_scattering(args.S02, args.energy_shift, 
                    args.absorber, args.ignore_elements, args.edge, args.rmax, 
                    trajectory)
        elif args.first_shell:
            k, chi = exafs_first_shell(args.S02, args.energy_shift, 
                    args.absorber, args.ignore_elements, args.edge, 
                    args.neighbor_cutoff, trajectory)
    
        return k, chi
