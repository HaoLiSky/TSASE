#!/usr/bin/env python
#lanczos algorithm for min-mode searching

import copy
import numpy as np
import sys
import os

from math import sqrt, atan, cos, sin, tan, pi
from ase  import atoms, units, io
from tsase.neb.util import vunit, vmag, vrand, sPBC
from tsase.dimer.ssdimer import SSDimer_atoms

class lanczos_atoms(SSDimer_atoms):

    def minmodesearch(self):
        F0    = self.update_general_forces(self.R0)
        F1    = self.rotation_update()
        dphi    = 100

        ## only for atom degree of freedoms
        #u = self.N[:-3].flatten()
        ## include cell too
        u = self.N.flatten()
        beta = vmag(u)
        size = (self.natom + 3) * 3
        #size = self.natom * 3
        T = np.zeros((size,size))
        Q = np.zeros((size,size))
            
        #while dphi > self.phi_tol and iteration < self.rotationMax:
        for i in range(size):
            Q[:, i] = u/beta
            Hv      = -(F1 - F0) / self.dR
            #u       = Hv[:-3].flatten()
            u       = Hv.flatten()
            if i > 0:
                u   = u - beta * Q[:,i-1]
            alpha   = np.vdot(Q[:, i], u)
            u       = u - alpha * Q[:, i]
            # re-orthogonalize
            u       = u - np.dot(Q, np.dot(Q.T, u))

            T[i, i] = alpha
            if i > 0:
                 T[i-1, i] = beta
                 T[i, i-1] = beta

            beta    = vmag(u)

            newN    = np.reshape(u/beta, (-1, 3))
            #self.N[:-3] = newN
            self.N  = newN
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
                evAtom_old = copy.copy(evAtom)
                evAtom  = np.reshape(evEst, (-1, 3))
                c0      = ew
                ## check with pi/2
                dotp = np.vdot(evAtom, evAtom_old)
                if dotp > 1: dotp = 1
                elif dotp < -1: dotp = -1
                dphi    = np.arccos(dotp) 
                if dphi > np.pi/2.0: dphi = np.pi - dphi
            else:
                #evAtom  = self.N[:-3]
                #c0      = np.vdot(Hv[:-3], evAtom)
                evAtom  = self.N
                c0      = np.vdot(Hv, evAtom)
            #print "i, Curvature, dphi:", i, c0, dphi
            if dphi < self.phi_tol or i == size-1 or i >= self.rotationMax:
                #self.N[:-3] = evAtom
                self.N = evAtom
                break

        self.curvature = c0
        return F0
        
