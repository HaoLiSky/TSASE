#!/usr/bin/env python
import numpy
import os
import ase
from ase import units
import string

def write_socorro(filename, atoms):
    R =  atoms.get_scaled_positions()
    box = atoms.get_cell()
    box /= units.Bohr

    f = open(filename, 'w')
    f.write('comment\n')
    f.write('1.0\n')
    for i in range(3):
        f.write('%f %f %f\n' % (box[i,0], box[i,1], box[i,2]))
    f.write('lattice\n')
    f.write('%i\n' % len(atoms))
    for i, atom in enumerate(atoms):
        f.write('%s %f %f %f\n' % (atom.symbol, R[i,0], R[i,1], R[i,2]) )

    f.close()

def read_socorro(filename):
    f = open(filename)
    comment = f.readline()
    lattice_constant = float(f.readline())
    cell = numpy.zeros((3,3))
    for i in range(3):
        line = f.readline()
        fields = line.split()
        for j in range(3):
            cell[i,j] = lattice_constant * float(fields[j]) * units.Bohr

    cartesian = True
    format = f.readline().strip()
    if format == 'lattice':
        cartesian = False

    natoms = int(f.readline())
    positions = numpy.zeros( (natoms, 3) )
    symbols = []
    for i, line in enumerate(f):
        fields = line.split()
        symbol = ''
        for ch in fields[0]:
            if ch in string.digits:
                break
            else:
                symbol += ch
        symbols.append(symbol)
        for j in range(3):
            positions[i,j] = float(fields[j+1])
            if cartesian:
                positions[i,j] *= lattice_constant

    atoms = ase.Atoms(symbols=symbols, cell=cell, pbc=True)
    if cartesian:
        atoms.set_positions(positions) * units.Bohr
    else:
        atoms.set_scaled_positions(positions)

    return atoms
