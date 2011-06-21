import numpy
import ase
import tsase

a = ase.io.read("al.traj")
al = tsase.calc.al()
a.set_calculator(al)

u = a.get_potential_energy()
print u, "eV"

f = a.get_forces()
print numpy.linalg.norm(f), "eV/angstrom"


