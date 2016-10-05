
import numpy as np
import lj_

class lj:

    def __init__(self, epsilon = 1.0, sigma = 1.0, cutoff = 10.0):
        self.sigma = sigma
        self.epsilon = epsilon
        self.cutoff = cutoff
        self.u0 = 4.0 * self.epsilon
        self.atoms = None
        self.u = None
        self.f = None
        self.force_calls = 0

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
        if self.f is None or self.u is None or atoms is None:
            return True
        return False

    def set_atoms(self, atoms):
        pass

    def calculate(self):
        ra = self.atoms.positions.ravel()
        fa = ra * 0.0
        uRet = np.array([0], 'd')
        ax = self.atoms.cell[0][0]
        ay = self.atoms.cell[1][1]
        az = self.atoms.cell[2][2]
        lj_.force(ra, fa, uRet, ax, ay, az, self.u0, self.epsilon, self.sigma, self.cutoff)
        self.f = np.resize(fa, (len(self.atoms),3))
        self.u = uRet[0]
        self.force_calls += 1


