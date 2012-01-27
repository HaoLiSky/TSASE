
import ase
import math
import numpy

def graphene(a=1.42, vacuum=8.0):
    y = math.sin(2*math.pi*120/360)
    x = math.cos(2*math.pi*120/360)
    g = numpy.array([[0.5, 0.0, 0.0],
                     [1.5, 0.0, 0.0],
                     [0.0, y,   0.0],
                     [2.0, y,   0.0]])
    p = ase.Atoms(symbols='CCCC', pbc=[True, True, True])
    p.cell[0][0] = 3 * a
    p.cell[1][1] = 2 * y * a
    p.cell[2][2] = vacuum
    p.positions = g * a
    p.center()
    return p


