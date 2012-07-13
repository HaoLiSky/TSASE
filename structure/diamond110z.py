
import numpy
from ase import Atoms

diamond110z_direct = numpy.array([[0.00, 0.00, 0.00],
                                  [0.50, 0.25, 0.00],
                                  [0.50, 0.50, 0.50],
                                  [0.00, 0.75, 0.50]])

def diamond110z(symbol, a = 5.467):
    p = Atoms(symbols=symbol*4, pbc=[True, True, True])
    p.cell[0][0] = (0.5**0.5) * a
    p.cell[1][1] = a
    p.cell[2][2] = (0.5**0.5) * a
    p.positions = numpy.dot(diamond110z_direct, p.cell)
    return p


