
import os
import tempfile
import shutil
import numpy as np

class bopfox:

    def __init__(self, atomsbx="atoms.bx", bondsbx="bonds.bx", infoxbx="infox.bx", bopfox="bopfox"):
        self.atoms = None
        self.atomsbx = atomsbx
        self.bondsbx = bondsbx
        self.infoxbx = infoxbx
        self.bopfox = bopfox

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

    def _write_fox(self, a):
        fout = open("struc.bx", 'w')
        fout.write("StrucName = struc\n")
        fout.write("aLat = 1.0\n")
        fout.write("a1 = % .12f   % .12f   % .12f\n" % (a.cell[0][0], a.cell[0][1], a.cell[0][2]))
        fout.write("a2 = % .12f   % .12f   % .12f\n" % (a.cell[1][0], a.cell[1][1], a.cell[1][2]))
        fout.write("a3 = % .12f   % .12f   % .12f\n" % (a.cell[2][0], a.cell[2][1], a.cell[2][2]))
        fout.write("coord = cartesian\n")
        symbols = a.get_chemical_symbols()
        r = a.get_positions()
        for i in range(len(a)):
            fout.write("%s % .12f % .12f % .12f /% 3d\n" % (symbols[i], r[i][0], r[i][1], r[i][2], i))
        fout.close()

    def _read_fu(self, N):
        fu = open("struc.EnFo.bx", 'r')
        u = float(fu.readline())
        f = np.zeros((N, 3))
        for i in range(N):
            line = fu.readline().strip().split()
            f[i][0] = float(line[0])
            f[i][1] = float(line[1])
            f[i][2] = float(line[2])
        fu.close()
        return f, u

    def calculate(self):
        curdir = os.getcwd()
        bfdir = tempfile.mkdtemp()
        os.system("cp %s %s %s %s" % (self.atomsbx, self.bondsbx, self.infoxbx, bfdir)) 
        os.chdir(bfdir)
        self._write_fox(self.atoms)
        os.system(self.bopfox)
        self.f, self.u = self._read_fu(len(self.atoms))
        os.chdir(curdir)
        shutil.rmtree(bfdir)



