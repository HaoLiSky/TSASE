import numpy
from numpy import exp

class ZDP_5Gauss():
	def __init__(self):
		self.atoms = None
		self.u = None
		self.f = None

	def calculate(self):
		ra = self.atoms.positions
		forces = ra * 0.0
		x = ra[0][0]
		y = ra[0][1]
		
		self.u =  -1.5 * exp(-((x - 1) ** 2) - ((y - 1) ** 2))
		self.u += -1.5 * exp(-((x - 3) ** 2) - ((y - 1) ** 2))
		self.u += -1.5 * exp(-((x - 1) ** 2) - ((y - 3) ** 2))
		self.u += -1.5 * exp(-((x - 3) ** 2) - ((y - 3) ** 2))
		self.u += -1.0 * exp(-3.0 * (((x - 2) ** 2) + ((y - 2) ** 2)))
		
		forces[0][0]  = 3. * exp(-((x - 1) ** 2) - ((y - 1) ** 2)) * (x - 1)
		forces[0][0] += 3. * exp(-((x - 1) ** 2) - ((y - 3) ** 2)) * (x - 1)
		forces[0][0] += 3. * exp(-((x - 3) ** 2) - ((y - 1) ** 2)) * (x - 3)
		forces[0][0] += 3. * exp(-((x - 3) ** 2) - ((y - 3) ** 2)) * (x - 3)
		forces[0][0] += 6. * exp(-3.0 * (((x - 2) ** 2) + ((y - 2) ** 2))) * (x - 2)
		
		forces[0][1]  = 3. * exp(-((x - 1) ** 2) - ((y - 1) ** 2)) * (y - 1)
		forces[0][1] += 3. * exp(-((x - 3) ** 2) - ((y - 1) ** 2)) * (y - 1)
		forces[0][1] += 3. * exp(-((x - 1) ** 2) - ((y - 3) ** 2)) * (y - 3)
		forces[0][1] += 3. * exp(-((x - 3) ** 2) - ((y - 3) ** 2)) * (y - 3)
		forces[0][1] += 6. * exp(-3.0 * (((x - 2) ** 2) + ((y - 2) ** 2))) * (y - 2)

		self.f = forces
        
	def get_potential_energy(self, atoms=None, force_consistent=False):
		if self.calculation_required(atoms, "energy"):
			self.atoms = atoms.copy()
			self.calculate()
		return self.u
		 
	def get_forces(self, atoms):
		if self.calculation_required(atoms, "forces"):
			self.atoms = atoms.copy()
			self.calculate()	
		return self.f.copy()
		                
	def get_stress(self, atoms):
		raise NotImplementedError
        
	def calculation_required(self, atoms, quantities):
		if atoms != self.atoms or self.atoms == None:
			return True
		if self.f == None or self.u == None or atoms == None:
			return True
		return False

	def set_atoms(self, atoms):
		pass

