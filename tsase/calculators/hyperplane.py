
#from pylab import *  ### 
import tsase
#import ase
import numpy
#import random
### Define function that converts V to H ############

class hyperplane_potential:

    def __init__(self,pot,negeig=None):
        self.pot = pot
        self.atoms = None
        self.min = self.pot.get_positions()
        if negeig == None:
            self.negeig = numpy.zeros(numpy.shape(p.get_positions()))
        else:
            self.negeig = negeig
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
            self.atoms = atoms.copy()
            diff = self.atoms.get_positions() - self.min
            dot = numpy.vdot(self.negeig,diff)
            newpos = self.atoms.get_positions() - dot*self.negeig
            self.atoms.set_positions(newpos)
            return True
        if self.f == None or self.u == None or atoms == None:
            self.atoms = atoms.copy()
            diff = self.atoms.get_positions() - self.min
            dot = numpy.vdot(self.negeig,diff)
            newpos = self.atoms.get_positions() - dot*self.negeig
            self.atoms.set_positions(newpos)
            return True
        else:
            return False

    def set_atoms(self, atoms):
        pass

    def potforce(self):
        u  = self.u1d()
        f  = self.f1d()
        return u, f

    def u1d(self):
        self.pot.set_positions(self.atoms.get_positions())
        return self.pot.get_potential_energy()

    def f1d(self):
        self.pot.set_positions(self.atoms.get_positions())
        f = self.pot.get_forces()
        f -= numpy.vdot(f,self.negeig)*self.negeig 
        return f


