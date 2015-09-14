#!/usr/bin/env python
'''
Calculate the smallest curvature of the isosurface. 
In the previous version, we use the curvature of the 
isosurface in the 2D plane spanned by self.N and F0 to 
estimate the samllest curevature. 
If it is positive, the dimer moves in the regular way; 
If it is zero, turn off the parallel force;
If it is negative, reverse the parallel force;

P. Xiao and Q. Wu and G. Henkelman Basin constrained k-dimer method for saddle point finding, J. Chem. Phys. 141 164111 (2014)     
'''

import time
import copy
import numpy as np
import sys
import os

from math import sqrt, atan, cos, sin, tan, pi
from ase  import atoms, units, io
from tsase.neb.util import vunit, vmag, vrand, sPBC
from tsase.io import read_con
from tsase.dimer.ssdimer import SSDimer_atoms
import sys

class KSSDimer_atoms(SSDimer_atoms):

    def __init__(self, R0 = None, mode = None, maxStep = 0.2, dT = 0.1, dR = 0.001, 
                 phi_tol = 5, rotationMax = 4, ss = False, express=np.zeros((3,3)), 
                 beta = 5, releaseF=0.1, phi_iso = 5,
                 rotationOpt = 'cg', weight = 1, noZeroModes = True):
        """
        gamma1     =  2.0 / (np.exp((self.isocurv - 0.0) * beta) + 1.0) - 1.0
        gamma2     = -1.0 / (np.exp((self.isocurv - 0.0) * beta) + 1.0) + 1.0
        self.Ftrans = gamma2 * Fperp - gamma1 * Fparallel
        %self.isocurv is kappa%

        New Parameters:
        beta       - how fast the forces change near the kappa=0 hyperplane 
        releaseF   - when the magnitude of total force below this value, switch to the regular dimer
        phi_iso    - rotation converging tolerence, in degree, for kappa search
        """
        SSDimer_atoms.__init__(self, R0, mode = mode, maxStep = maxStep, dT = dT, dR = dR, 
                 phi_tol = phi_tol, rotationMax = rotationMax, ss = ss, express = express, 
                 rotationOpt = rotationOpt, weight = weight, noZeroModes = noZeroModes)
        self.alpha    = 0.0
        self.beta     = beta
        self.releaseF = releaseF
        self.phi_iso  = phi_iso /180.0 * pi
        self.Contour  = self.N.copy() #mode perpendicular to the force, along the contour(isosurface)
   

    def get_forces(self):
        F0 = self.minmodesearch(minoriso = 'min')
        Fparallel = np.vdot(F0, self.N) * self.N
        Fperp     = F0 - Fparallel
        beta1      = self.beta
        if self.curvature > 0:
            self.Ftrans = -Fparallel
            print "drag up directly"
        else:
            if vmag(F0) < self.releaseF: 
                gamma1 = 1; gamma2 = 1; 
            else:
                self.minmodesearch(minoriso = 'iso')
                print "exact isocurve:", self.isocurv
                if abs(self.isocurv) > 20: self.isocurv = 20 * np.sign(self.isocurv)
                gamma1     =  2.0 / (np.exp((self.isocurv - 0.0) * beta1) + 1.0) - 1.0
                gamma2     = -1.0 / (np.exp((self.isocurv - 0.0) * beta1) + 1.0) + 1.0
            self.Ftrans = gamma2 * Fperp - gamma1 * Fparallel

            #####xph: needs further improvement to insure k-dimer won't be trapped at isocurve=0.0 and Fperp=0.0
            #if vmag(F0) > 0.03 and vmag(self.Ftrans) < 0.005:
                #self.Ftrans = Fperp - Fparallel
                #self.Ftrans *= 10

        if self.ss:
            return self.Ftrans
        else:
            return self.Ftrans[:-3]
 
    def rotation_update(self, mode):
        # position of R1
        self.iset_endpoint_pos(mode, self.R0, self.R1)
        # force of R1
        F1    = self.update_general_forces(self.R1)
        return F1

    def rotation_plane(self, Fperp, Fperp_old, mode):
        # determine self.T (big theta in the paper) with CG method
        a = abs(np.vdot(Fperp, Fperp_old))
        b = np.vdot(Fperp_old, Fperp_old)
        if a <= 0.5*b and b != 0:  
            gamma = np.vdot(Fperp, Fperp-Fperp_old) / b
        else:
            gamma = 0
        Ttmp       = Fperp + gamma * self.T * self.Tnorm
        Ttmp       = Ttmp - np.vdot(Ttmp, mode) * mode
        self.Tnorm = np.linalg.norm(Ttmp)
        self.T     = vunit(Ttmp)

    def minmodesearch(self, minoriso):
        # rotate dimer to the minimum mode direction
        # self.N, the dimer direction; 
        # self.T, the rotation direction, spans the rotation plane with self.N.

        ## self.N or self.Contour
        if minoriso == "min":
            self.F0 = self.update_general_forces(self.R0)
            mode    = self.N.copy()
        else:
            F0unit  = vunit(self.F0)
            mode    = self.Contour.copy()
            mode    = mode - np.vdot(mode, F0unit) * F0unit

        # project out any rigid translation and rotation
        if not self.ss: mode    = self.project_translt_rott(mode, self.R0)

        F0    = self.F0
        F1    = self.rotation_update(mode)

        phi_min = 1.5
        Fperp   = F1 * 0.0 # avoid Fperp_old asignment error
        iteration = 0
            
        while iteration < self.rotationMax:
            
            ## self.N or self.Contour
            if minoriso == "min":
                if phi_min <= self.phi_tol: break
                F0perp    = F0 - np.vdot(F0, mode) * mode
                F1perp    = F1 - np.vdot(F1, mode) * mode
            else:
                if phi_min <= self.phi_iso: break
                mode      = mode - np.vdot(mode, F0unit) * F0unit
                F0perp    = F0 * 0.0
                F1iso     = F1 - np.vdot(F1, F0unit) * F0unit
                F1perp    = F1iso - np.vdot(F1iso, mode) * mode

            Fperp_old = Fperp
            Fperp     = 2.0 * (F1perp - F0perp)
            if iteration == 0: 
                Fperp_old = Fperp
                self.T    = mode * 0.0
                self.Tnorm= 0.0
            # update self.T
            self.rotation_plane(Fperp, Fperp_old, mode)
            if minoriso == "iso": 
                self.T = self.T - np.vdot(self.T, F0unit) * F0unit
                self.T = vunit(self.T)

            # project out any rigid translation and rotation
            if not self.ss: self.T = self.project_translt_rott(self.T, self.R0)

            # curvature and its derivative
            c0     = np.vdot(F0-F1, mode) / self.dR
            c0d    = np.vdot(F0-F1, self.T) / self.dR * 2.0
            phi_1  = -0.5 * atan(c0d / (2.0 * abs(c0)))

            if abs(phi_1) <= self.phi_tol: break
            # calculate F_prime: force after rotating the dimer by phi_prime
            N1_prime = vunit(mode * cos(phi_1) + self.T * sin(phi_1))
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
            mode = vunit(mode * cos(phi_min) + self.T * sin(phi_min))
            c0 = c0_min
            # update F1 by linear extropolation
            F1 = F1 * (sin(phi_1 - phi_min) / sin(phi_1)) + F1_prime * (sin(phi_min) / sin(phi_1)) \
                 + F0 * (1.0 - cos(phi_min) - sin(phi_min) * tan(phi_1 * 0.5))
            #print "phi_min:", phi_min, "toque:", vmag(Fperp)
            iteration += 1

        if minoriso == "min":
            self.curvature = c0
            self.N         = mode.copy()
            '''
            ## estimate the curvature of the isophote (isoenergy curve)
            Fg = vunit(F0)
            cosalpha = abs(np.vdot(Fg, self.N)) 
            print "cosalpha", cosalpha
            if cosalpha > 0.01 and cosalpha < 0.999:
                N1_prime = vunit(self.N - np.vdot(Fg, self.N)*Fg)
                self.iset_endpoint_pos(N1_prime, self.R0, self.R1_prime)
                F1_prime = self.update_general_forces(self.R1_prime)
                c0_prime = np.vdot(F0-F1_prime, N1_prime) / self.dR 
                self.kappa = -c0_prime/vmag(F0)
            else: 
                self.kappa = -0.25
            '''
        else:
            self.isocurv   = -c0/vmag(F0)
            self.Contour   = mode.copy()
            print "isocurv iteration:", iteration
        return F0
        

