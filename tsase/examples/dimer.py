import ase
from ase import dimer
import tsase

atoms = tsase.io.read_con("al.con")
al = tsase.calculators.al()
atoms.set_calculator(al)
mma = dimer.MinModeAtoms(atoms)
opt = dimer.MinModeTranslate(mma, trajectory = "dimer.traj")
opt.run(fmax=0.01)

