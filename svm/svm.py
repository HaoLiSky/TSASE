import numpy
from numpy import exp

try:
	import sklearn
except:
	import scikits.learn

class SVM_dynamics():
	def __init__(self,calc,selSVM,k):
		
		self.atoms = None
		self.calc = calc
		self.selSVM = selSVM
		self.u = None
		self.f = None
		self.k = k
	
	def set_k(self,k):
		self.k = k
		
	def get_k(self):
		return self.k

	def calculate(self):
		self.calc.get_forces(self.atoms)
		self.f = self.calc.f.copy()
		datapoint = numpy.asarray([self.atoms.get_positions().ravel()])
		grad = self.svm_rbf_deriv(self.selSVM, datapoint) 
		gradnorm = numpy.sqrt(numpy.vdot(grad,grad))
		try:
			grad /= gradnorm # the gradient is now normalized
			prediction = self.selSVM.decision_function(datapoint)[0,0] / gradnorm
		except:
			prediction = 0.
		self.f += -1 * self.k * prediction * numpy.reshape(grad,numpy.shape(self.atoms.get_positions()))
		
	def svm_rbf_deriv(self,clf,Points): # Don't forget this returns the positive of the gradient!
		import numpy as np
		try:
			from sklearn.metrics.pairwise import rbf_kernel
		except:
			from scikits.learn.metrics.pairwise import rbf_kernel
		Grad = np.zeros(Points.shape)
		noP = Points.shape[0]
		alpha = clf.dual_coef_
		d = np.ones((len(clf.support_),1))
		# loop over test points an compute the gradient individually
		for k in range(noP):
			P = Points[k,:]
			P.shape = 1,-1
			# be careful with the gamma paremeter: gamma= 1/(2sigma^2)
			if (float(sklearn.__version__) <= 0.8):
				Kmatrix = rbf_kernel(clf.support_vectors_,P,np.sqrt(1./(clf.gamma*2)))
			else: 
				Kmatrix = rbf_kernel(clf.support_vectors_,P,clf.gamma)
			kstar = alpha.T *Kmatrix	
			b = (clf.support_vectors_-np.dot(d,P))*2*clf.gamma
			Grad[k,:] = np.dot(kstar.T,b)
		Grad *= -1 # The function returns the negative of the gradient without this change
		return Grad
		
	def get_potential_energy(self, atoms, force_consistent=False):
		if self.calculation_required(atoms, "energy"):
			self.atoms = atoms.copy()
			self.calc.atoms = atoms.copy()
			self.calc.calculate()
		return self.calc.u
		 
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

