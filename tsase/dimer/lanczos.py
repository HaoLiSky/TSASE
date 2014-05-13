#!/usr/bin/env python

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

class lanczos_atoms:

    def __init__(self, R0 = None, mode = None, maxStep = 0.2, dT = 0.1, dR = 0.005, 
                 phi_tol = 10, rotationMax = 4, ss = True, express=np.zeros((3,3)), 
                 weight = 1, lowestmode = None):
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
        self.ss       = ss
        self.express  = express
   
        vol           = self.R0.get_volume()
        avglen        = (vol/self.natom)**(1.0/3.0)
        self.weight   = weight
        self.jacobian = avglen * self.natom**0.5 * self.weight
        self.V = np.zeros((self.natom+3,3))
        # for comparison
        self.lowestmode = lowestmode
        if lowestmode == None:
            print 'please provide the true lowest mode for comparison'
            sys.exit()
        self.record = [[], [], []]

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
        F0 = self.minmodesearch()
        Fparallel = np.vdot(F0, self.N) * self.N
        #if self.curvature > 0 and self.steps < 50:
        if self.curvature > 0:
            self.Ftrans = -Fparallel
            print "drag up directly"
        else:
            self.Ftrans = F0 - 2.0 * Fparallel
        return self.Ftrans
 
    def iset_endpoint_pos(self, Ni, R0, Ri, dRi):
        # update the position of Ri
        dRvec = dRi * Ni
        #cell0 = R0.get_cell()
        #cell1 = cell0 + np.dot(cell0, dRvec[-3:]) / self.jacobian
        #Ri.set_cell(cell1, scale_atoms=True)
        #vdir  = R0.get_scaled_positions()
        #ratom = np.dot(vdir, cell1) + dRvec[:-3]
        ratom = R0.get_positions() + dRvec[:-3]
        Ri.set_positions(ratom)
 
    def rotation_update(self):
        # position of R1
        self.iset_endpoint_pos(self.N, self.R0, self.R1, self.dR)
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
        Ttmp       = Ttmp - np.vdot(Ttmp, self.N) * self.N
        self.Tnorm = np.linalg.norm(Ttmp)
        self.T     = vunit(Ttmp)
        
    def minmodesearch(self):
        # rotate dimer to the minimum mode direction
        # self.N, the dimer direction; 
        # self.T, the rotation direction, spans the rotation plane with self.N.
        F0    = self.update_general_forces(self.R0)
        F1    = self.rotation_update()
        dphi    = 100

        #r = self.N
        ##only for atom degree of freedoms
        u = self.N[:-3].flatten()
        beta = vmag(u)
        #size = (self.natom + 3) * 3
        size = self.natom * 3
        T = np.zeros((size,size))
        Q = np.zeros((size,size))
            
        #while dphi > self.phi_tol and iteration < self.rotationMax:
        for i in range(size):
            Q[:, i] = u/beta
            Hv      = -(F1 - F0) / self.dR
            u       = Hv[:-3].flatten()
            if i > 0:
                u   = u - beta * Q[:,i-1]
            alpha   = np.vdot(Q[:, i], u)
            u       = u - alpha * Q[:, i]
            u       = u - np.dot(Q, np.dot(Q.T, u))

            T[i, i] = alpha
            if i > 0:
                 T[i-1, i] = beta
                 T[i, i-1] = beta

            beta    = vmag(u)

            newN    = np.reshape(u/beta, (-1, 3))
            ## check with pi/2
            dphi    = np.arccos(np.vdot(newN, self.N[:-3])) 
            if dphi > np.pi/2.0: dphi = np.pi - dphi
            self.N[:-3] = newN
            F1      = self.rotation_update()

            ##Check Eigenvalues
            if i > 0:
                eigenValues, eigenVectors  = np.linalg.eig(T[:i+1, :i+1])
                idx = eigenValues.argsort()   
                ew  = eigenValues[idx][0]
                evT = eigenVectors[:,idx][:, 0]
                ##Convert eigenvector of T matrix to eigenvector of full Hessian
                evEst = Q[:, :i+1].dot(evT)
                evEst = vunit(evEst)
                evAtom  = np.reshape(evEst, (-1, 3))
                c0      = ew
                theta2true = np.arccos(np.vdot(evAtom, self.lowestmode[:-3]))
                if theta2true > np.pi/2.0: theta2true = np.pi - theta2true
                self.record[0].append(i)
                self.record[1].append(theta2true)
                '''
                ##check orthonality
                nonorth = 0
                for j in range(i):
                    nonorth += abs(np.vdot(Q[:, i], Q[:, j]))
                nonorth /= i
                print nonorth
                #self.record[2].append(nonorth)
                '''
                self.record[2].append(c0)
                ##only for comparison
                if theta2true < self.phi_tol: 
                    break

            else:
                evAtom  = self.N[:-3]
                c0      = np.vdot(Hv[:-3], evAtom)
            #print "i, Curvature:", i, c0
            '''
            if dphi < self.phi_tol or i == size-1:
                self.N[:-3] = evAtom
                break
            '''

        self.curvature = c0
        print "Lanczos nstpes:", self.forceCalls
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
            Rcenter = self.R0.copy()
            Rshadow = self.R0.copy()
            self.iset_endpoint_pos(self.N, Rcenter, Rshadow, 0.5)
            Oshadow = Rshadow[0] #copy the oxygen atom 
            Rcenter.append(Oshadow)
            io.write(movie, Rcenter, format='vasp')
        # While the max atom force is greater than some criteria...
        while self.getMaxAtomForce() > minForce and self.forceCalls < maxForceCalls:
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
                                   
