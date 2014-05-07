
import numpy
import morse_

class morse():
    def __init__(self, u0 = 0.7102, alpha = 1.6047, r0 = 2.8970, rc = 9.5):
        self.u0 = u0        # well depth
        self.alpha = alpha  # stiffness of potential
        self.r0 = r0        # diatomic distance
        self.rc = rc        # cutoff distance
        self.atoms = None
        self.u = None
        self.f = None
        self.force_calls = 0

    def calculate(self):
        ra = self.atoms.positions.ravel()
        fa = ra * 0.0
        uRet = numpy.array([0.0])
        ax = self.atoms.cell[0][0]
        ay = self.atoms.cell[1][1]
        az = self.atoms.cell[2][2]
        morse_.force(ra, fa, uRet, ax, ay, az, self.u0, self.alpha, self.r0, self.rc)
        self.f = numpy.resize(fa, (len(self.atoms),3))
        self.u = uRet[0]
        self.force_calls += 1
        
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

        
