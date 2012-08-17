
import ase
from tsase.data import elements
import math
import numpy

def randomCluster(n):
    atoms = ase.Atoms()
    atoms.append(ase.Atom('H', [0,0,0]))
    for i in range(n-1):
        print i
        farthest = 1e300
        a = numpy.random.normal(0, 1, 3)
        a = n * a/numpy.linalg.norm(a)
        A = numpy.array([0,0,0]) - a
        Au = A/numpy.linalg.norm(A)
        A = 2*n*Au
        for atom in atoms:
            b = atom.position
            B = b - a
            C = Au * numpy.dot(B, Au)
            D = b - (a + C)
            Dm = numpy.linalg.norm(D)
            Q = elements['H']['radius'] + elements[atom.symbol]['radius']
            if Dm > Q:
                continue
            if Dm == 0:
                V = Q
            else:            
                V = math.sqrt((Q*Q) - (Dm*Dm))
            farthest = min(numpy.linalg.norm(C)-V, farthest)
        p = a + Au * farthest                    
        atoms.append(ase.Atom('H', p))
    return atoms
        
    
    
