
import numpy
from numpy import cos, sin, pi

class voter97():
	def __init__(self, d1=4.0, d2=1.0, d3=0., d4=1.0):
		self.d1=d1
		self.d2=d2
		self.d3=d3
		self.d4=d4
		self.k=2.0 * pi
		self.atoms = None
		self.u = None
		self.f = None

	def calculate(self):
		ra = self.atoms.positions
		forces = ra * 0.0
		a = 1.0+self.d1*(ra[0][1]-1.0)
		kx = self.k*ra[0][0]
		ky = self.k*(ra[0][1]-1.0)
		forces[0][0] = self.k*a*sin(kx)+self.d3*self.k*sin(self.k*ra[0][0]/self.d4)/self.d4
		forces[0][1] = -self.d2*self.k*ky-self.d1*cos(kx)
		self.f = forces
		self.u = cos(kx) * a + self.d2/2* (ky ** 2) + self.d3*cos(kx/self.d4)
        
	def get_potential_energy(self, atoms, force_consistent=False):
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

