#!/usr/bin/env python
import numpy
import os
from tsase.io import write_socorro

class Socorro:
    def __init__(self):
        self.forces = None
        self.energy = None
        self.force_calls = 0
        self.atoms = None

    def set_atoms(self, atoms):
        pass

    def calculation_required(self, atoms, quantities):
        if self.energy is None and self.atoms is not None:
            return True
        if atoms != self.atoms:
            return True
        if (atoms.get_positions()!=self.atoms.get_positions()).any():
            return True
        return False

    def get_potential_energy(self, atoms=None, force_consistent=False):
        if self.calculation_required(atoms, 'energy'):
            self.atoms = atoms.copy()
            self.calculate_forces(atoms)
        return self.energy

    def get_forces(self, atoms):
        if self.calculation_required(atoms, 'forces'):
            self.atoms = atoms.copy()
            self.calculate_forces(atoms)
        return self.forces

    def get_stress(self, atoms):
        raise NotImplementedError

    def run(self):
        pass

    def calculate_forces(self, atoms):
        write_socorro('crystal', atoms)

        cmd = os.environ['SOCORRO_COMMAND']
        os.system('%s > %s' % (cmd, 'stdout'))

        self.force_calls += 1
        energy, forces = parse_diary('diaryf')

        self.energy = energy
        self.forces = forces

def parse_diary(filename):
    f = open(filename)

    force_section = False
    force_section_count = 0
    forces = []
    for line in f:
        line = line.strip()
        fields = line.split()
        if line.startswith('cell energy'):
            energy = float(fields[3])
        if line.startswith('Atomic forces:'): 
            force_section = True

        if force_section:
            force_section_count += 1
            if force_section_count <= 3: continue
            if len(line) == 0:
                force_section = False
                continue
            forces.append([ float(fields[1]), float(fields[2]), float(fields[3]) ])

    forces = numpy.array(forces)
    return energy, forces
