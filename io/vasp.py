import ase
import numpy

def read_xdatcar(fileName):
    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()
    lattice_constant = float(lines[1].strip())
    cell = numpy.array([[float(x) * lattice_constant for x in lines[2].split()], 
                        [float(x) * lattice_constant for x in lines[3].split()], 
                        [float(x) * lattice_constant for x in lines[4].split()]])
    elements = lines[5].split()
    natoms = [int(x) for x in lines[6].split()]
    nframes = (len(lines)-7)/(sum(natoms) + 1)
    trajectory = []
    for i in range(nframes):
        a = ase.Atoms()
        a.cell = cell.copy()
        j = 0
        for N, e in zip(natoms, elements):
            for k in range(N):
                x = float(lines[8 + i * (sum(natoms) + 1) + j].split()[0])
                y = float(lines[8 + i * (sum(natoms) + 1) + j].split()[1])
                z = float(lines[8 + i * (sum(natoms) + 1) + j].split()[2])
                x, y, z = cell.dot(numpy.array([x,y,z]))
                a.append(ase.Atom(e, (x, y, z), mass = 1.0))
                j += 1
        trajectory.append(a)
    return trajectory
