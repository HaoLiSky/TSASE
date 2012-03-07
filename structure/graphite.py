
import ase
import math
import numpy

def graphite(a=1.42, interlayer=3.35):
    y = math.sin(2*math.pi*120/360)
    g = numpy.array([[0.5, 0.0, 0.0],
                     [1.5, 0.0, 0.0],
                     [0.0, y,   0.0],
                     [2.0, y,   0.0],
                     [1.5, 0.0, 1.0],
                     [2.5, 0.0, 1.0],
                     [1.0, y,   1.0],
                     [3.0, y,   1.0]])
    p = ase.Atoms(symbols='CCCCCCCC', pbc=[True, True, True])
    p.cell[0][0] = 3 * a
    p.cell[1][1] = 2 * y * a
    p.cell[2][2] = 2 * interlayer
    p.positions[:,0:2] = g[:,0:2] * a
    p.positions[:,2] = g[:,2] * interlayer
    p.center()
    return p


