import ase
import al

pot = al.eam_al_calculator()
p = ase.io.read("POSCAR")
p.set_calculator(pot)
opt = ase.optimize.MDMin(p, trajectory = "al.traj")
opt.run(fmax = 0.0001)

