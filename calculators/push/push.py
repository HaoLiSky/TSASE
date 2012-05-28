
import numpy as np

class push:

    def __init__(self, attractive_force=0.1, repulsive_force=1.0, bond_distance=2.5):
        self.attractive_force = attractive_force
        self.repulsive_force = repulsive_force
        self.bond_distance = bond_distance
        self.atoms = None

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

    def calculate(self):
        rf = self.repulsive_force
        af = self.attractive_force
        z = self.bond_distance
        r2 = rf/2
        r22 = (rf**2)/2
        r24 = (rf**2)/4
        a2 = af/2
        a22 = (af**2)/2
        a24 = (af**2)/4
        f = np.zeros((len(self.atoms), 3))
        u = 0
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                if i == j:
                    continue
                # Get the vector and distance between the two atoms.
                v = self.atoms.positions[i] - self.atoms.positions[j]
                vr = np.linalg.solve(self.atoms._cell.T, v)
                v = np.dot(vr - np.round(vr) * self.atoms._pbc, self.atoms._cell)
                d = np.linalg.norm(v)  
                if d == 0:
                    raise ValueError("Push potential cannot handle zero distance between atoms.")
                # The repulsive regime:
                if d < -r2 + z:
                    u += -rf*(d-z) - r22 + r24
                    f[i] += rf * (v/d)
                # The attractive regime:
                elif d > a2 + z:
                    u += af*(d-z) - a22 + a24
                    f[i] += -af * (v/d)
                # The harmonic regime:
                else:
                    u += (d-z)**2
                    f[i] += (-2*(d-z)) * (v/d)
        self.f = f
        self.u = u / 2.0



