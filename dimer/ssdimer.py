#!/usr/bin/env python

import time
import copy
import numpy as np
import sys
import os

from math import sqrt, atan, cos, sin, tan, pi
from ase  import atoms, units, io
from tsase.neb.util import vunit, vmag, vrand, sPBC

class SSDimer_atoms:

    def __init__(self, R0 = None, mode = None, maxStep = 0.2, dT = 0.1, dR = 0.005, 
                 phi_tol = 10, rotationMax = 4, express=np.zeros((3,3)), weight = 1):
        """
        Parameters:
        force - the force to use
        R0 - starting point
        mode - initial mode (will be randomized if one is not provided)
        maxStep - longest distance dimer can move in a single iteration
        dT - quickmin timestep
        dR - finite difference step size for translation
        phi_tol - rotation converging tolerence, degree
        rotationMax - max rotations per translational step
        """
        self.steps = 0
        self.dT = dT
        self.dR = dR
        self.phi_tol = phi_tol /180.0 * pi
        self.R0    = R0
        self.N = mode
        self.natom  = self.R0.get_number_of_atoms()
        if self.N == None:
            self.N = vrand(np.zeros((self.natom+3,3)))
        self.N = vunit(self.N)
        self.maxStep = maxStep
        self.Ftrans = None
        self.forceCalls = 0
        self.R1 = self.R0.copy()
        self.R1_prime = self.R0.copy()
        calc    = self.R0.get_calculator()
        self.R1.set_calculator(calc)
        self.R1_prime.set_calculator(calc)
        self.rotationMax = rotationMax
        self.express  = express
   
        vol           = self.R0.get_volume()
        avglen        = (vol/self.natom)**(1.0/3.0)
        self.weight   = weight
        self.jacobian = avglen * self.natom**0.5 * self.weight
        self.V = np.zeros((self.natom+3,3))

    # Pipe all the stuff from Atoms that is not overwritten.
    # Pipe all requests for get_original_* to self.atoms0.
    def __getattr__(self, attr):
        """Return any value of the Atoms object"""
        return getattr(self.R0, attr)
        print "*************************************"
        print attr

    def __len__(self):
        return self.natom+3

    def get_positions(self):
        r    = self.R0.get_positions()*0.0
        Rc   = np.vstack((r, self.R0.get_cell()*0.0))
        return Rc

    def set_positions(self,dr):
        rcell  = self.R0.get_cell()
        rcell += np.dot(rcell, dr[-3:]) / self.jacobian
        self.R0.set_cell(rcell, scale_atoms=True)
        ratom  = self.R0.get_positions() + dr[:-3]
        self.R0.set_positions(ratom)
    
    def update_general_forces(self, Ri):
        #update the generalized forces (f, st)
        self.forceCalls += 1
        f    = Ri.get_forces()
        stt  = Ri.get_stress()
        vol  = Ri.get_volume()*(-1)
        st   = np.zeros((3,3))
        #following the order of get_stress in vasp.py
        #(the order of stress in ase are the same for all calculators)
        st[0][0] = stt[0] * vol  
        st[1][1] = stt[1] * vol
        st[2][2] = stt[2] * vol
        st[2][1] = stt[3] * vol
        st[2][0] = stt[4] * vol
        st[1][0] = stt[5] * vol
        st  -= self.express * (-1)*vol
        #print "original stress (no projecton applied):"
        #print st
        Fc   = np.vstack((f, st/self.jacobian))
        return Fc
  
    def get_curvature(self):
        return self.curvature

    def get_forces(self):
        F0 = self.minmodesearch()
        Fparallel = np.vdot(F0, self.N) * self.N
        if self.curvature > 0:
            self.Ftrans = -Fparallel
        else:
            self.Ftrans = F0 - 2.0 * Fparallel
        return self.Ftrans
 
    def iset_endpoint_pos(self, Ni, R0, Ri):
        # update the position of Ri
        dRvec = self.dR * Ni
        cell0 = R0.get_cell()
        cell1 = cell0 + np.dot(cell0, dRvec[-3:]) / self.jacobian
        Ri.set_cell(cell1, scale_atoms=True)
        ratom = R0.get_positions() + dRvec[:-3]
        Ri.set_positions(ratom)
 
    def rotation_update(self):
        # position of R1
        self.iset_endpoint_pos(self.N, self.R0, self.R1)
        # force of R1
        F1    = self.update_general_forces(self.R1)
        return F1
        
    def rotation_plane(self, Fperp, Fperp_old):
        # determine self.T (big theta in the paper) with CG method
        a = abs(np.vdot(Fperp, Fperp_old))
        b = np.vdot(Fperp_old, Fperp_old)
        if a <= 0.5*b and b != 0:  
            gamma = np.vdot(Fperp, Fperp-Fperp_old) / b
        else:
            gamma = 0
        Ttmp       = Fperp + gamma * self.T * self.Tnorm
        self.Tnorm = np.linalg.norm(Ttmp)
        self.T     = vunit(Ttmp)
        
    def minmodesearch(self):
        # rotate dimer to the minimum mode direction
        # self.N, the dimer direction; 
        # self.T, the rotation direction, spans the rotation plane with self.N.


        F0    = self.update_general_forces(self.R0)
        F1    = self.rotation_update()
        F0perp = F0 - np.vdot(F0, self.N) * self.N

        phi_min = 1.5
        Fperp   = F1 * 0.0 # avoid Fperp_old asignment error
        iteration = 0
        while abs(phi_min) > self.phi_tol and iteration < self.rotationMax:

            F1perp    = F1 - np.vdot(F1, self.N) * self.N
            Fperp_old = Fperp
            Fperp     = 2.0 * (F1perp - F0perp)
            if iteration == 0: 
                Fperp_old = Fperp
                self.T    = self.N * 0.0
                self.Tnorm= 0.0
            # update self.T
            self.rotation_plane(Fperp, Fperp_old)
            # curvature and its derivative
            c0     = np.vdot(F0-F1, self.N) / self.dR
            c0d    = np.vdot(F0-F1, self.T) / self.dR * 2.0
            phi_1  = -0.5 * atan(c0d / (2.0 * abs(c0)))

            if abs(phi_1) <= self.phi_tol: break
            # calculate F_prime: force after rotating the dimer by phi_prime
            N1_prime = vunit(self.N * cos(phi_1) + self.T * sin(phi_1))
            self.iset_endpoint_pos(N1_prime, self.R0, self.R1_prime)
            F1_prime = self.update_general_forces(self.R1_prime)
            c0_prime = np.vdot(F0-F1_prime, N1_prime) / self.dR 
            
            # calculate phi_min
            b1 = 0.5 * c0d
            a1 = (c0 - c0_prime + b1 * sin(2 * phi_1)) / (1 - cos(2 * phi_1))
            a0 = 2 * (c0 - a1)
            phi_min = 0.5 * atan(b1 / a1)
            c0_min  = 0.5 * a0 + a1 * cos(2.0 * phi_min) + b1 * sin(2 * phi_min)

            # check whether it is minimum or maximum
            if c0_min > c0 :
                phi_min += pi * 0.5
                c0_min   = 0.5 * a0 + a1 * cos(2.0 * phi_min) + b1 * sin(2 * phi_min)
                 
            # update self.N
            self.N = vunit(self.N * cos(phi_min) + self.T * sin(phi_min))
            # update F1 by linear extropolation
            F1 = F1 * (sin(phi_1 - phi_min) / sin(phi_1)) + F1_prime * (sin(phi_min) / sin(phi_1)) \
                 + F0 * (1.0 - cos(phi_min) - sin(phi_min) * tan(phi_1 * 0.5))
            iteration += 1
        self.curvature = c0
        return F0
        

#######################################################################################################
# The following part can be replaced by FIRE or MDMin optimizer in ase, see the ssdimer.py in examples
    def step(self):
        self.steps += 1
        Ftrans = self.get_forces() 
        print "Ftrans",vmag(Ftrans)
        dV = Ftrans * self.dT
        if np.vdot(self.V, Ftrans) > 0:
            self.V = dV * (1.0 + np.vdot(dV, self.V) / np.vdot(dV, dV))
        else:
            self.V = dV
        step = self.V * self.dT
        if vmag(step) > self.maxStep:
            step = self.maxStep * vunit(step)
 
        self.set_positions(step)
        self.E = self.get_potential_energy()
    
    def getMaxAtomForce(self):
        if self.Ftrans is None:
            return 1000
        maxForce = -1
        #for i in range(len(self.FCross)):
        #    maxForce = max(maxForce, vmag(self.FCross[i]))
        maxForce = vmag(self.Ftrans)
        return maxForce
            
    def search(self, minForce = 0.01, quiet = False, maxForceCalls = 100000, movie = None, interval = 50):
        self.converged = False
        if movie:
            #savexyz(movie, self.R0, 'w')
            io.write(movie, self.R0, format='vasp')
        # While the max atom force is greater than some criteria...
        while self.getMaxAtomForce() > minForce and self.forceCalls < maxForceCalls:
            # Take a Dimer step.
            self.step()
            if movie and self.steps % interval == 0:
                io.write('movie.tmp', self.R0, format='vasp')
                os.system('cat movie.tmp >> '+movie)
            if not quiet:
                ii = self.steps
                ff = self.getMaxAtomForce()
                cc = self.curvature
                ee = self.E
                nf = self.forceCalls
                if self.steps % 100 == 0 or self.steps == 1:
                    print "Iteration       Force       Curvature        Energy     ForceCalls"
                    print "-------------------------------------------------------------------------------"
                    print "%3i %13.6f %13.6f %13.6f %3i" % (ii,float(ff),float(cc),float(ee),nf)
                else:
                    print "%3i %13.6f %13.6f %13.6f %3i" % (ii,float(ff),float(cc),float(ee),nf)
        if self.getMaxAtomForce() <= minForce:
            self.converged = True
                                   
