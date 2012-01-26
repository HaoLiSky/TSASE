from ase.constraints import FixBondLengths

def FixBlock(indices_list):
    bonds = []
    for indices in indices_list:
        for i in range(len(indices) - 1):
            for j in range(i + 1, len(indices)):
                bonds.append([indices[i], indices[j]])
    return FixBondLengths(bonds)
    



    