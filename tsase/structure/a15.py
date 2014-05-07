
import numpy
from ase import Atoms

a15_direct = numpy.array([[0.00, 0.00, 0.00],
                          [0.50, 0.50, 0.50],
                          [0.25, 0.50, 0.00],
                          [0.75, 0.50, 0.00],
                          [0.00, 0.25, 0.50],
                          [0.00, 0.75, 0.50],
                          [0.50, 0.00, 0.25],
                          [0.50, 0.00, 0.75]])

def a15(symbol, a = 5.0295756):
    if a == None:
        a = 5.0295756
    p = Atoms(symbols=symbol*8, pbc=[True, True, True])
    p.cell[0][0] = a
    p.cell[1][1] = a
    p.cell[2][2] = a
    p.positions = a * a15_direct
    return p


