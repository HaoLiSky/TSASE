from ase.constraints import FixBondLengths

def FixBlock(indices):
    bonds = []
    for i in range(len(indices) - 1):
        for j in range(i + 1, len(indices)):
            bonds.append([i, j])
    return FixBondLengths(bonds)
    



    