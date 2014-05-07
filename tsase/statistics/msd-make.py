#!/usr/bin/env python
import ase
from ase.io import read
from expectra.io import read_xdatcar
import numpy as np
import sys
import math

def main():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--ignore-elements', type=str, metavar='ELEMENTS',
            help='comma delimited list of elements to ignore in the ' + \
            'pdf')
    parser.add_argument('--skip', type=int, default=0,
            help='number of frames to skip at the beginning')
    parser.add_argument('--every', type=int, default=1,
            help='number of frames to between each step')
    parser.add_argument('--timestep', type=float, default=2,
            help='time-step used in MD simulation')
    parser.add_argument('trajectories', metavar='TRAJ', nargs='+',
            help='trajectory file (POSCAR, con, xyz)')

    args = parser.parse_args()

    if args.ignore_elements:
        args.ignore_elements = args.ignore_elements.split(',')

    snapshots_command(args)

def drange(start, stop, step):
    while start < stop:
        yield start
        start += step

def snapshots_command(args):
    trajectory = []
    for filename in args.trajectories:
        print 'reading', filename
        if filename[-3:] == 'con':
            trajectory += read(filename, type=eon)[args.skip::args.every]
        elif filename.startswith('POSCAR') or filename.startswith('CONTCAR'):
            trajectory = [read(filename)]
        else:
            trajectory += read_xdatcar(filename, args.skip, args.every)
#    print trajectory
    df= args.timestep*float(args.every)
    tmax=len(trajectory)*df
    cnt = np.zeros(len(trajectory))
    ldt = np.zeros(len(trajectory))
    atoms = trajectory[0]
    natoms = len(atoms)
    sdx = np.zeros(shape=(len(trajectory)))
    sdy = np.zeros(shape=(len(trajectory)))
    sdz = np.zeros(shape=(len(trajectory)))
    msd = np.zeros(shape=(len(trajectory)))
    #print len(sdx[:,0])   #print tmax
    #print f
    ## Loop over snapshots
    for i, atoms in enumerate(trajectory):
        atoms = atoms.copy()
#        cell =  atoms.get_cell()
#        bl = cell[2,2]
#        bl_norm = np.linalg.norm(cell)
#        bin_size= bl_norm/2.0/(args.num_bin)
        if args.ignore_elements:
           ignore_indicies = [atom.index for atom in atoms 
                   if atom.symbol in args.ignore_elements]
           del atoms[ignore_indicies]
        pos_t = atoms.get_positions()

    ###Loop over snapshots

        for j, atoms2 in enumerate(trajectory):
            if j > i:
                atoms2 = atoms2.copy()
                if args.ignore_elements:
                    ignore_indicies = [atom.index for atom in atoms   
                        if atom.symbol in args.ignore_elements]
                    del atoms[ignore_indicies]
                pos_dt = atoms2.get_positions()
                dt = (j-i)
         #       print dt
                for k in xrange(len(atoms)):
                        cnt[dt] += 1
                        sdx[dt] += (pos_dt[k,0] - pos_t[k,0])*(pos_dt[k,0] - pos_t[k,0])
                        sdy[dt] += (pos_dt[k,1] - pos_t[k,1])*(pos_dt[k,1] - pos_t[k,1])
                        sdz[dt] += (pos_dt[k,2] - pos_t[k,2])*(pos_dt[k,2] - pos_t[k,2])

    #print cnt
    for d in xrange(1, len(trajectory)):
     #   print sdx[d]
        sdx[d] /= cnt[d]/natoms
        sdy[d] /= cnt[d]/natoms
        sdz[d] /= cnt[d]/natoms
        msd[d] = np.sqrt(sdx[d]+sdy[d]+sdz[d])
        ldt[d] = float(d)*float(df)
    print ldt
    ## WRITE TO FILE ##
    outdat=np.column_stack((ldt, sdx, sdy, sdz, msd))
    np.savetxt("msd.dat", outdat, delimiter="   ")

 #   outdat_interp=np.column_stack((xnew, ynew, yint))
 #   np.savetxt("interprdf.dat", outdat_interp, delimiter="   ")


if __name__ == '__main__':
        main()

