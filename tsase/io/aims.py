import numpy
import ase

def read_aims(filename):
    f = open(filename)

    atoms = ase.Atoms()

    l = 0
    cell = numpy.zeros( (3,3) )
    for line in f:
        line = line[:line.find('#')]
        line = line.strip()
        if len(line) == 0: continue
        fields = line.split()

        if fields[0] == 'lattice_vector':
            cell[l,0] = float(fields[1])
            cell[l,1] = float(fields[2])
            cell[l,2] = float(fields[3])
            l += 1
            if l == 3:
                atoms.set_cell(cell)
                print cell

        if fields[0] == 'atom':
            pos = numpy.array((float(fields[1]),float(fields[2]),float(fields[3])))
            element = fields[4]
            atoms.append(ase.Atom(element, position=pos))

        if fields[0] == 'atom_frac':
            scaled_pos = numpy.array((float(fields[1]),float(fields[2]),float(fields[3])))
            pos = numpy.dot(cell, scaled_pos)
            element = fields[4]
            atoms.append(ase.Atom(element, position=pos))
    f.close()
    return atoms

def write_aims(filename, atoms):
    f = open(filename, 'w')

    for atom in atoms:
        f.write("atom %.10f %.10f %.10f %s\n" % (atom.position[0], 
            atom.position[1],atom.position[2], atom.symbol) )
    f.close()
