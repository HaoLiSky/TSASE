import numpy as np
import tsase

class BGSD:
    def __init__(self,pot,alpha=1.0,beta=1.0,dr=1e-3):
        self.alpha = alpha
        self.beta = beta
        self.pot = pot
        self.dr = dr
        self.atoms = None
        self.u = None
        self.f = None

    def calculate(self):
        self.u, self.f = self.potforce()

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

    def potforce(self):      
        u  = self.u1d()
        f  = self.f1d()
        return u, f
        
    def u1d(self):
        self.pot.set_positions(self.atoms.get_positions())
        f = self.pot.get_forces()
        g2 = np.vdot(f,f)
        u = self.pot.get_potential_energy()
        har = self.alpha * (u - self.beta) * (u - self.beta)
        toteng = g2 + har
        return toteng
 
    def f1d(self):
        self.dr = 0.001
        self.pot.set_positions(self.atoms.get_positions())
        f = self.pot.get_forces()
        u = self.pot.get_potential_energy()
        pos = self.pot.get_positions()
        gradhar = 2*self.alpha*(u - self.beta)*f
        self.pot.set_positions(pos + self.dr*f)
        gradf = - (self.pot.get_forces() - f)/self.dr
        self.pot.set_positions(pos)
        totf = gradf +gradhar
        return totf
