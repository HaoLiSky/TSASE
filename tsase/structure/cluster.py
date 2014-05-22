
import ase
from tsase.data import elements
import math
import random
import numpy

def randomCluster(n):
    atoms = ase.Atoms()
    atoms.append(ase.Atom('H', [0,0,0]))
    for i in range(n-1):
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
        
    
def fastRandomCluster(n, gridSize=1.0, rattle=0.1):

    def getNeighborSlots(position):
        options = []
        if (position[0] + gridSize, position[1] + 0, position[2] + 0) not in taken:
            options.append((position[0] + gridSize, position[1] + 0, position[2] + 0))
        if (position[0] - gridSize, position[1] + 0, position[2] + 0) not in taken:
            options.append((position[0] - gridSize, position[1] + 0, position[2] + 0))
        if (position[0] + 0, position[1] + gridSize, position[2] + 0) not in taken:
            options.append((position[0] + 0, position[1] + gridSize, position[2] + 0))
        if (position[0] + 0, position[1] - gridSize, position[2] + 0) not in taken:
            options.append((position[0] + 0, position[1] - gridSize, position[2] + 0))
        if (position[0] + 0, position[1] + 0, position[2] + gridSize) not in taken:
            options.append((position[0] + 0, position[1] + 0, position[2] + gridSize))
        if (position[0] + 0, position[1] + 0, position[2] - gridSize) not in taken:
            options.append((position[0] + 0, position[1] + 0, position[2] - gridSize))
        return options

    atoms = ase.Atoms()
    taken = {}
    available = [(0,0,0)]
    for i in range(n):
        selected = random.choice(available)
        taken[selected] = True
        available.remove(selected)
        newSlots = getNeighborSlots(selected)
        for slot in newSlots:
            if slot not in available:
                available.append(slot)
        atoms.append(ase.Atom('H', selected))
    dr = numpy.random.normal(0, 1, (n,3))
    big = max([numpy.linalg.norm(a) for a in dr])
    for a in dr:
        a /= big
    dr *= rattle
    atoms.positions += dr
    return atoms


if __name__ == "__main__":
    from ase.optimize import FIRE
    import ase.io
    import tsase
    p = fastRandomCluster(100)
    p.center(1000)
    lj = tsase.calculators.lj()
    p.set_calculator(lj)
    dyn = FIRE(p, trajectory='relax.traj')
    dyn.run(fmax=0.05)


