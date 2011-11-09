import tsase

p = tsase.io.read_bopfox("wa15.bx")

bf = tsase.calculators.bopfox("atoms.bx", "bonds.bx", "infox.bx", "bopfox")
p.set_calculator(bf)

print p.get_forces()
print p.get_potential_energy()
