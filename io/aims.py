import numpy
import ase

def read_aims(filename):
    f = open(filename)

    atoms = ase.Atoms()

    for line in f:
        line = line[:line.find('#')]
        line = line.strip()
        if len(line) == 0: continue
        fields = line.split()

        if fields[0] == 'atom':
            pos = numpy.array((float(fields[1]),float(fields[2]),float(fields[3])))
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
