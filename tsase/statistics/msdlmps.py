#!/usr/bin/env python
import ase
from ase.io import read
from expectra.io import read_xdatcar
from tsase.io.lammps import read_dump
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
            help='number of frames between each step')
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
        elif filename.startswith('trj'):
            trajectory += read_dump(filename)
        else:
            trajectory += read_xdatcar(filename, args.skip, args.every)

    ### Set variables and array sizes
    df= args.timestep*args.every
    tmax=len(trajectory)*df
    cnt = np.zeros(len(trajectory))
    ldt = np.zeros(len(trajectory))
    atoms = trajectory[0]
    sdx = np.zeros(shape=(len(trajectory)))
    sdy = np.zeros(shape=(len(trajectory)))
    sdz = np.zeros(shape=(len(trajectory)))
    msd = np.zeros(shape=(len(trajectory)))
    natoms = 0
    ###Loop over snapshots-outer
    for i, atoms in enumerate(trajectory):
    #    if i == 0:
        atoms = atoms.copy()

        if args.ignore_elements:
           ignore_indicies = [atom.index for atom in atoms 
                   if atom.symbol in args.ignore_elements]
           del atoms[ignore_indicies]
        pos_sc = atoms.get_scaled_positions()
        cell= atoms.get_cell()
        pos_t = np.dot(pos_sc, cell)

    ###Loop over snapshots-inner--track dt
        for j, atoms2 in enumerate(trajectory):
            if j > i:
                atoms2 = atoms2.copy()

                if args.ignore_elements:
                    ignore_indicies = [atom.index for atom in atoms
                        if atom.symbol in args.ignore_elements]
                    del atoms[ignore_indicies]
                pos_sc = atoms2.get_scaled_positions()
                
                pos_dt = np.dot(pos_sc, atoms2.get_cell())
                dt = (j-i)
               
                cnt[dt] += 1
            ### get difference squared
                for k in xrange(len(atoms)):
                        if j == 1:
                            natoms += 1
                        sdx[dt] += (pos_dt[k,0] - pos_t[k,0])*(pos_dt[k,0] - pos_t[k,0])
                        sdy[dt] += (pos_dt[k,1] - pos_t[k,1])*(pos_dt[k,1] - pos_t[k,1])
                        sdz[dt] += (pos_dt[k,2] - pos_t[k,2])*(pos_dt[k,2] - pos_t[k,2])

    ## normalize to # dt diff and # atoms
    print cnt
    print natoms
    for d in xrange(1, len(trajectory)):
        sdx[d] /= (cnt[d]*natoms)
        sdy[d] /= (cnt[d]*natoms)
        sdz[d] /= (cnt[d]*natoms)
        msd[d] = np.sqrt(sdx[d]+sdy[d]+sdz[d])
        ldt[d] = float(d)*float(df)

    ## WRITE TO FILE ##
    outdat=np.column_stack((ldt, sdx, sdy, sdz, msd))
    np.savetxt("msd.dat", outdat, delimiter="   ")


if __name__ == '__main__':
        main()

