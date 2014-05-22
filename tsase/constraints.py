from ase.constraints import FixBondLengths
import numpy as np

def FixBlock(indices_list):
    bonds = []
    for indices in indices_list:
        for i in range(len(indices) - 1):
            for j in range(i + 1, len(indices)):
                bonds.append([indices[i], indices[j]])
    return FixBondLengths(bonds)

class FixBlockAxes:

    def __init__(self, index):
        self.index = index

    def adjust_positions(self, oldpositions, newpositions):
        step = newpositions[self.index] - oldpositions[self.index]
        step = step.sum(0)/len(step)
        newpositions[self.index] = oldpositions[self.index] + np.array([step])

    def adjust_forces(self, positions, forces):
        forces[self.index] = forces[self.index]*0 + forces[self.index].sum(0)/len(self.index)

    def copy(self):
        return FixBlockAxes(self.index)


