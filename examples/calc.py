import numpy
import ase
import tsase

a = ase.io.read("al.traj")
al = tsase.calculators.al()
a.set_calculator(al)

u = a.get_potential_energy()
print u, "eV"

f = a.get_forces()
print numpy.linalg.norm(f), "eV/angstrom"

b=[a]
b+=[a.copy() for i in range(2)]
r = a.get_positions()
r[0][0]+=1
b[1].set_positions(r)

u = b[0].get_potential_energy()
print u, "eV"
b[1].set_calculator(b[0].get_calculator())
u = b[1].get_potential_energy()
print u, "eV"
u = b[2].get_potential_energy()
print u, "eV"
