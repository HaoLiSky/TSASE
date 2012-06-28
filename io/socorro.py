#!/usr/bin/env python
import numpy
import os
import ase
from ase import units

def write_socorro(filename, atoms):
    N = atoms.get_number_of_atoms()
    Z = atoms.get_atomic_numbers()

    R =  atoms.get_positions()
    R /= units.Rydberg
    box = atoms.get_cell()
    box /= units.Rydberg

    f = open(filename, 'w')
    f.write('comment\n')
    f.write('1.0\n')
    for i in range(3):
        f.write('%f %f %f\n' % (box[i,0], box[0,1], box[0,2]))
    f.write('cartesian\n')
    f.write('%i\n' % len(atoms))
    for atom in atoms:
        f.write('%s %f %f %f\n' % ( atom.symbol, atom.position[0], 
                atom.position[1], atom.position[2]))
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
        symbols.append(fields[0])
        for j in range(3):
            positions[i,j] = float(fields[j+1]) * units.Bohr
            if cartesian:
                positions[i,j] *= lattice_constant

    atoms = ase.Atoms(symbols=symbols, cell=cell, pbc=True)
    if cartesian:
        atoms.set_positions(positions)
    else:
        atoms.set_scaled_positions(positions)

    return atoms
