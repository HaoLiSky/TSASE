#!/usr/bin/env python
'''
Calculate the smallest curvature of the isosurface. 
In the previous version, we use the curvature of the 
isosurface in the 2D plane spanned by self.N and F0 to 
estimate the samllest curevature. 
If it is positive, the dimer moves in the regular way; 
If it is zero, turn off the parallel force;
If it is negative, reverse the parallel force;
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
from tsase import neb
import sys

class SSDimer_atoms:

    def __init__(self, R0 = None, mode = None, maxStep = 0.2, dT = 0.1, dR = 0.005, 
                 phi_tol = 10, phi_iso = 10, rotationMax = 4, ss = True, express=np.zeros((3,3)), 
                 estimateF1 = True, nebInitiate = False, originalRotation = False, 
                 dTheta = 1, beta = 5, releaseF=0.1, weight = 1, noZeroModes = True):
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
        ss - boolean, solid-state dimer or regular dimer. Default: ssdimer
        """
        self.steps = 0
        self.dT = dT
        self.dR = dR
        self.phi_tol = phi_tol / 180.0 * pi
        self.phi_iso = phi_iso / 180.0 * pi
        self.R0    = R0
        self.N = mode
        self.natom  = self.R0.get_number_of_atoms()
        if self.N == None:
            self.N = vrand(np.zeros((self.natom+3,3)))
        self.N = vunit(self.N)
        self.Contour = self.N.copy() #mode perpendicular to the force, along the contour(isosurface)
        self.maxStep = maxStep
        self.Ftrans = None
        self.forceCalls = 0
        self.R1 = self.R0.copy()
        self.R1_prime = self.R0.copy()
        calc    = self.R0.get_calculator()
        self.R1.set_calculator(calc)
        self.R1_prime.set_calculator(calc)
        self.rotationMax = rotationMax
        self.noZeroModes = noZeroModes # Set to False for 2D model potentials
        self.ss       = ss
        self.express  = express
        self.estimateF1  = estimateF1
        self.nebInitiate = nebInitiate
        self.originalRotation = originalRotation
        self.dTheta   = dTheta / 180.0 * pi
        self.alpha    = 0.0
        self.beta     = beta
        self.releaseF = releaseF
   
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
        r    = self.R0.get_positions()
        if self.ss:
            # return zeros, so the vector passed to set_positions is just dr in the generalized space. 
            # otherwise, "+" and "-" operations for position vectors need to be redefined in the optimizer
            # this trick only works for first order optimzers for sure, where no operation applied to position vectors until the last update.
            Rc   = np.vstack((r*0.0, self.R0.get_cell()*0.0))
            return Rc
        else:
            return r

    def set_positions(self,dr):
        if self.ss:
            rcell  = self.R0.get_cell()
            rcell += np.dot(rcell, dr[-3:]) / self.jacobian
            self.R0.set_cell(rcell, scale_atoms=True)
            ratom  = self.R0.get_positions() + dr[:-3]
            self.R0.set_positions(ratom)
        else:
            # get_positions() returns non-zero values
            # thus this dr is the final positions, not just dr
            ratom  = dr
            self.R0.set_positions(ratom)
    
    def update_general_forces(self, Ri):
        #update the generalized forces (f, st)
        self.forceCalls += 1
        f    = Ri.get_forces()
        if self.ss: stt  = Ri.get_stress()
        vol  = Ri.get_volume()*(-1)
        st   = np.zeros((3,3))
        #following the order of get_stress in vasp.py
        #(the order of stress in ase are the same for all calculators)
        if self.ss:
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

    def get_mode(self):
        return self.N

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
                if abs(self.isocurv) > 20: self.isocurv = 20 * np.sign(self.isocurv)
                gamma1     =  2.0 / (np.exp((self.isocurv - 0.0) * beta1) + 1.0) - 1.0
                gamma2     = -1.0 / (np.exp((self.isocurv - 0.0) * beta1) + 1.0) + 1.0
                #gamma1     =  1.2 / (np.exp((self.isocurv - 0.0) * beta1 + np.log(5)) + 1.0) - 0.2
                #gamma1     =  1.5 / (np.exp((self.isocurv - 0.0) * beta1 ) + 1.0) - 0.5
                #gamma2     = (-1.0 / (np.exp((self.isocurv + 0.0) * beta1) + 1.0) + 1.0) * 0.5
            self.Ftrans = gamma2 * Fperp - gamma1 * Fparallel

            #####xph: needs further improvement to insure k-dimer won't be trapped at isocurve=0.0 and Fperp=0.0
            #if vmag(F0) > 0.02 and vmag(self.Ftrans) < 0.01 and self.alpha == 0.0:
            #    self.alpha = 1.0
            #if self.alpha == 1.0:
            #    self.Ftrans = Fperp - Fparallel
            #####or:
            #if vmag(F0) > 0.03 and vmag(self.Ftrans) < 0.005:
                #self.Ftrans = Fperp - Fparallel
                #self.Ftrans *= 10

            print "exact isocurve:", self.isocurv
        return self.Ftrans
 
    def project_translt_rott(self, N, R0):
        if not self.noZeroModes: return N
        # Project out rigid translational mode
        for axisx in range(3):
            transVec = np.zeros((self.natom+3, 3))
            transVec[:-3, axisx] = 1.0
            transVec = vunit(transVec)
            N -= np.vdot(N, transVec)*transVec
        # Project out rigid rotational mode
        for axisx in ['x', 'y', 'z']:
            ptmp = R0.copy()
            # rotate a small angle around the center of mass
            ptmp.rotate(axisx, 0.02, center='COM', rotate_cell=False)
            rottVec = ptmp.get_positions() - R0.get_positions()
            rottVec = vunit(rottVec)
            #if np.vdot(N[:-3], rottVec) > 0.1: print "mainly rotation around "+axisx
            N[:-3] -= np.vdot(N[:-3], rottVec)*rottVec
        return N
 
    def iset_endpoint_pos(self, Ni, R0, Ri, dRi):
        # update the position of Ri
        dRvec = dRi * Ni
        cell0 = R0.get_cell()
        cell1 = cell0 + np.dot(cell0, dRvec[-3:]) / self.jacobian
        Ri.set_cell(cell1, scale_atoms=True)
        vdir  = R0.get_scaled_positions()
        ratom = np.dot(vdir, cell1) + dRvec[:-3]
        ratom = R0.get_positions() + dRvec[:-3]
        Ri.set_positions(ratom)
 
    def rotation_update(self, mode):
        # position of R1
        self.iset_endpoint_pos(mode, self.R0, self.R1, self.dR)
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
            self.iset_endpoint_pos(N1_prime, self.R0, self.R1_prime, self.dR)
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
            if self.estimateF1:
                F1 = F1 * (sin(phi_1 - phi_min) / sin(phi_1)) + F1_prime * (sin(phi_min) / sin(phi_1)) \
                     + F0 * (1.0 - cos(phi_min) - sin(phi_min) * tan(phi_1 * 0.5))
            else: 
                F1    = self.rotation_update(mode)
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
                self.iset_endpoint_pos(N1_prime, self.R0, self.R1_prime, self.dR)
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
        

#######################################################################################################
# The following part can be replaced by FIRE or MDMin optimizer in ase, see the ssdimer.py in examples
    def step(self):
        self.steps += 1
        Ftrans = self.get_forces() 
        #print "Ftrans",vmag(Ftrans)
        dV = Ftrans * self.dT
        #if np.vdot(self.V, Ftrans) > 0 and self.isocurv <= 0:
        if np.vdot(self.V, Ftrans) > 0 :
            self.V = dV * (1.0 + np.vdot(dV, self.V) / np.vdot(dV, dV))
        else:
            self.V = dV
        step = self.V * self.dT
        if vmag(step) > self.maxStep:
            step = self.maxStep * vunit(step)
 
        self.set_positions(self.get_positions() + step)
        self.E = self.get_potential_energy()
        
    
    def getMaxAtomForce(self):
        if self.Ftrans is None:
            return 1000
        maxForce = -1
        #for i in range(len(self.FCross)):
        #    maxForce = max(maxForce, vmag(self.FCross[i]))
        maxForce = vmag(self.Ftrans)
        return maxForce
            
    def search(self, minForce = 0.01, quiet = False, maxForceCalls = 300000, movie = None, interval = 50):
        self.converged = False
        if movie:
            #savexyz(movie, self.R0, 'w')
            Rcenter = self.R0.copy()
            Rshadow = self.R0.copy()
            self.iset_endpoint_pos(self.N, Rcenter, Rshadow, 0.5)
            Oshadow = Rshadow[0] #copy the oxygen atom 
            Rcenter.append(Oshadow)
            io.write(movie, Rcenter, format='vasp')
        # While the max atom force is greater than some criteria...
        try: cc = self.curvature
        except: cc = 1.0
        while (self.getMaxAtomForce() > minForce or cc > 0) and self.forceCalls < maxForceCalls :
            # Take a Dimer step.
            self.step()

            if movie and self.steps % interval == 0:
                Rcenter = self.R0.copy()
                Rshadow = self.R0.copy()
                self.iset_endpoint_pos(self.N, Rcenter, Rshadow, 0.5)
                Oshadow = Rshadow[0] #copy the oxygen atom 
                Rcenter.append(Oshadow)
                io.write('movie.tmp', Rcenter, format='vasp')
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
                                   
