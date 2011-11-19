"""
The parallel version of generalized nudged elastic path (ssneb) module.
Multi processing scheme is employed here. When evaluating the force, each process deals with one image.
For mpirun and vasp, the execution in run_vasp.py should be like:
mpirun -np 8 --hostfile my_hosts vasp
Only works for SGE job system now, because my_hosts is greped from  $PE_HOSTFILE 
The number of hosts used should be equal to the number of images.
"""

import numpy
import os,sys
from copy import deepcopy
from math import sqrt, atan, pi
from util import vmag, vunit, vproj, vdot, sPBC
from ase import atoms
from multiprocessing import Pool

def calpot((imgi,p,host)):
        # making a directory for each image, which is nessary for vasp to read last step's WAVECAR  
        # also, it is good to prevent overwriting files for parallelizaiton over images
        fdname = str('0'+str(imgi))
        if not os.path.exists(fdname): os.mkdir(fdname)
        os.chdir(fdname)
        hostfile = open('my_hosts','w')
        #hostfile.write(self.hosts[imgi]+' slots=8 \n')
        hostfile.write(host+' \n')
        hostfile.close()
        p.u = p.get_potential_energy()
        p.f = p.get_forces()
        stt = p.get_stress()
        os.chdir("../")
        p.icell = numpy.linalg.inv(p.get_cell())
        p.vdir  = p.get_scaled_positions()
        try: 
            p.st
        except:
            p.st = numpy.zeros((3,3))
        vol = p.get_volume()*(-1)
        p.st[0][0] = stt[0] * vol
        p.st[1][1] = stt[1] * vol
        p.st[2][2] = stt[2] * vol
        p.st[2][1] = stt[3] * vol
        p.st[2][0] = stt[4] * vol
        p.st[1][0] = stt[5] * vol
        p.st[0][1] = 0.0
        p.st[0][2] = 0.0
        p.st[1][2] = 0.0
        return p
     

class pssneb:
    """
    The generalized nudged elastic path (ssneb) class.
    """
    def __init__(self, p1, p2, numImages = 7, k = 5.0, tangent = "new",       \
                 dneb = False, dnebOrg = False, method = 'normal',            \
                 onlyci = False, weight = 1, parallel = False):
        """
        The neb constructor.
        Parameters:
            p1.......... one endpoint of the path
            p2.......... the other endpoint of the path
            numImages... the total number of images in the path, including the 
                         endpoints
            k........... the spring force constant
            tangent..... the tangent method to use, "new" for the new tangent,
                         anything else for the old tangent
            dneb........ set to true to use the double-nudging method
            dnebOrg..... set to true to use the original double-nudging method
            method...... "ci" for the climbing image method, anything else for
                         normal NEB method 
        """
        
        self.numImages = numImages
        self.k = k * numImages
        self.tangent = tangent
        self.dneb = dneb
        self.dnebOrg = dnebOrg
        self.method = method
        self.onlyci = onlyci
        self.weight = weight 
        self.parallel = parallel 

        #check the orientation of the cell, make sure a is along x, b is on xoy plane
        for p in [p1,p2]:
            cr = p.get_cell()
            if cr[0][1]**2+cr[0][2]**2+cr[1][2]**2 > 1e-3: 
                print "check the orientation of the cell, make sure a is along x, b is on xoy plane"
                sys.exit()
                
        #set the path by linear interpolation between end points
        n = self.numImages - 1
        self.path = [p1]
        self.path+= [p1.copy() for i in range(self.numImages-2)]
        self.path+= [p2]
        cell1 = p1.get_cell()
        cell2 = p2.get_cell()
        dRB   = (cell2 - cell1) / n # path for cell
        #don't use get_scaled_positions() or applay sPBC() here 
        #if the atoms can move over half of the lattice from initail to final
        icell = numpy.linalg.inv(cell1)
        vdir1 = numpy.dot(p1.get_positions(),icell)
        icell = numpy.linalg.inv(cell2)
        vdir2 = numpy.dot(p2.get_positions(),icell)
        dR    = sPBC(vdir2 - vdir1) / n # path for direct coordinates
        calc  = p1.get_calculator()
        for i in range(1, n):
            cellt = cell1 + dRB * i
            vdirt = vdir1 + dR * i
            rt    = numpy.dot(vdirt,cellt)
            self.path[i].set_cell(cellt)
            self.path[i].set_positions(rt)
            self.path[i].set_calculator(calc)
        self.Umaxi = 1

        #calculat the Jacobian to make cell move has the same unit and weight as atom move
        vol1     = self.path[0].get_volume()
        vol2     = self.path[self.numImages-1].get_volume()
        vol      = (vol1+vol2)*0.5
        self.natom = len(self.path[0]) 
        avglen   = (vol/self.natom)**(1.0/3.0)
        self.jacobian = avglen * self.natom**0.5 * self.weight
       
        #parse hostfile
        os.system("cat $PE_HOSTFILE > tmp")
        hfile = open('tmp','r')
        lines = hfile.readlines()
        self.hosts = []
        for line in lines:
            self.hosts.append(line.split()[0]) 
        if len(self.hosts) < self.numImages -2:
            print 'Warning: number of hosts < number of middle images'
            raw_input('If you are not using mpirun, hit Enter to continue')

        #add some new properties
        pool = Pool(processes=2)              # start several worker processes
        img0 = (0, self.path[0], self.hosts[0])
        img1 = (n, self.path[n], self.hosts[1])
        self.path[0], self.path[n] = pool.map(calpot, [img0, img1])         
        pool.close()
        self.path[0].cellt  = self.path[0].get_cell() * self.jacobian
        self.path[n].cellt  = self.path[n].get_cell() * self.jacobian

    def forces(self):
        """
        Calculate the forces for each image on the path.  Applies the force due
        to the potential and the spring forces.
        Parameters:
            force - the potential energy force.
        """
        
        # Calculate the force due to the potential on the intermediate points.
        pool = Pool(processes=self.numImages-2)              # start several worker processes
        # writing input and do the calculation in images' directories respectly
        mid_images = zip(range(1,self.numImages-1), self.path[1:-1], self.hosts)
        self.path[1:-1] = pool.map(calpot, mid_images)         
        pool.close()

        self.Umax  = self.path[0].u
        self.Umaxi = 0
        for i in range(1, self.numImages - 1):
            if self.path[i].u > self.Umax:
                self.Umax  = self.path[i].u
                self.Umaxi = i
            
        for i in range(1, self.numImages - 1):
            self.path[i].cellt  = self.path[i].get_cell() * self.jacobian
        # Loop over each intermediate point and calculate the tangents.
        for i in range(1, self.numImages - 1):

            # Here st should be cauchy stress tensor times cell volume. 
            # Timing box volume should have been done.
            self.path[i].totalf = numpy.vstack((self.path[i].f, self.path[i].st / self.jacobian))
            # realtf that needed by nebspline.pl is saved for output 
            self.path[i].realtf = self.path[i].totalf
            
            # If we're using the 'old' tangent, the tangent is defined as the
            # vector from the point behind the current image to the point in
            # front of the current image.
            # Haven't implemented for ssneb
            if self.tangent == 'old':
                self.path[i].n = (self.path[i + 1].r - self.path[i - 1].r)
            
            # Otherwise, we're using the 'new' tangent.
            # Ref:
            # G. Henkelman and H. Jonsson,  Improved tangent estimate in the 
            # nudged elastic path method for finding minimum energy paths and 
            # saddle points, J. Chem. Phys. 113, 9978-9985 (2000)
            else:
            
                # it wouldn't hurt for these names to be more descriptive
                UPm1 = self.path[i - 1].u > self.path[i].u
                UPp1 = self.path[i + 1].u > self.path[i].u
                
                # if V(i+1)>V(i)>V(i-1)
                # or V(i+1)<V(i)<V(i-1)
                # (this is the usual along the MEP)
                '''
                tangent
                '''
                if(UPm1 != UPp1):
                    if(UPm1):
                        # use direct coordinates to avoid double counting cell's movement
                        dr_dir  = sPBC(self.path[i].vdir - self.path[i - 1].vdir)
                        avgbox  = 0.5*(self.path[i].get_cell() + self.path[i - 1].get_cell())
                        sn  = numpy.dot(dr_dir,avgbox)
                        dh  = self.path[i].cellt - self.path[i - 1].cellt
                        snb = numpy.dot(self.path[i].icell, dh)*0.5 + numpy.dot(self.path[i - 1].icell, dh)*0.5
                        #---------------another way to average strain----------------------
                        #iavgbox = numpy.linalg.inv(avgbox)
                        #snb = numpy.dot(iavgbox, snb)
                        #------------------------------------------------------------------
                        self.path[i].n = numpy.vstack((sn,snb))
                    else:
                        dr_dir  = sPBC(self.path[i + 1].vdir - self.path[i].vdir)
                        avgbox  = 0.5*(self.path[i+1].get_cell() + self.path[i].get_cell())
                        sn  = numpy.dot(dr_dir,avgbox)
                        dh  = self.path[i + 1].cellt - self.path[i].cellt
                        snb = numpy.dot(self.path[i].icell, dh)*0.5 + numpy.dot(self.path[i + 1].icell, dh)*0.5
                        #---------------another way to average strain----------------------
                        #iavgbox = numpy.linalg.inv(avgbox)
                        #snb = numpy.dot(iavgbox, snb)
                        #------------------------------------------------------------------
                        self.path[i].n = numpy.vstack((sn,snb))
                # otherwise, we are near some extremum 
                else:
                    Um1 = self.path[i - 1].u - self.path[i].u
                    Up1 = self.path[i + 1].u - self.path[i].u
                    Umin = min(abs(Up1), abs(Um1))
                    Umax = max(abs(Up1), abs(Um1))
                    if(Um1 > Up1):
                        dr_dir  = sPBC(self.path[i + 1].vdir - self.path[i].vdir)
                        avgbox  = 0.5*(self.path[i + 1].get_cell() + self.path[i].get_cell())
                        sn      = numpy.dot(dr_dir,avgbox) * Umin
                        dr_dir  = sPBC(self.path[i].vdir - self.path[i - 1].vdir)
                        avgbox  = 0.5*(self.path[i].get_cell() + self.path[i - 1].get_cell())
                        sn     += numpy.dot(dr_dir,avgbox) * Umax

                        dh   = self.path[i + 1].cellt - self.path[i].cellt
                        snb1 = numpy.dot(self.path[i].icell, dh)*0.5 + numpy.dot(self.path[i + 1].icell, dh)*0.5
                        dh   = self.path[i].cellt - self.path[i - 1].cellt
                        snb2 = numpy.dot(self.path[i].icell, dh)*0.5 + numpy.dot(self.path[i - 1].icell, dh)*0.5
                        snb  = snb1 * Umin + snb2 * Umax
                        self.path[i].n = numpy.vstack((sn,snb))
                    else:
                        dr_dir  = sPBC(self.path[i + 1].vdir - self.path[i].vdir)
                        avgbox  = 0.5*(self.path[i + 1].get_cell() + self.path[i].get_cell())
                        sn      = numpy.dot(dr_dir,avgbox) * Umax
                        dr_dir  = sPBC(self.path[i].vdir - self.path[i - 1].vdir)
                        avgbox  = 0.5*(self.path[i].get_cell() + self.path[i - 1].get_cell())
                        sn     += numpy.dot(dr_dir,avgbox) * Umin

                        dh   = self.path[i + 1].cellt - self.path[i].cellt
                        snb1 = numpy.dot(self.path[i].icell, dh)*0.5 + numpy.dot(self.path[i + 1].icell, dh)*0.5
                        dh   = self.path[i].cellt - self.path[i - 1].cellt
                        snb2 = numpy.dot(self.path[i].icell, dh)*0.5 + numpy.dot(self.path[i - 1].icell, dh)*0.5
                        snb  = snb1 * Umax + snb2 * Umin
                        self.path[i].n = numpy.vstack((sn,snb))
        
        # Normalize the tangents.
        for i in range(1,self.numImages-1):
            self.path[i].n = vunit(self.path[i].n)

        # Loop over each intermediate image and adjust the potential energy 
        # force and apply the spring force.
        for i in range(1, self.numImages - 1):
        
            # Push the climbing image uphill.
            if self.method == 'ci' and i == self.Umaxi:
                self.path[i].totalf -= 2.0 * vproj(self.path[i].totalf, self.path[i].n) 
                self.path[i].fPerp   = self.path[i].totalf
            
            # And for the non-climbing images...
            else:
                
                # Calculate the force perpendicular to the tangent. 
                self.path[i].fPerp = self.path[i].totalf - vproj(self.path[i].totalf,   \
                                                            self.path[i].n)
                # Calculate the spring force.
                Rm1  = sPBC(self.path[i - 1].vdir - self.path[i].vdir)
                avgbox  = 0.5*(self.path[i - 1].get_cell() + self.path[i].get_cell())
                Rm1  = numpy.dot(Rm1,avgbox) 
                dh   = self.path[i - 1].cellt - self.path[i].cellt
                Rm1b = numpy.dot(self.path[i].icell, dh)*0.5 + numpy.dot(self.path[i - 1].icell, dh)*0.5
                Rm1  = numpy.vstack((Rm1,Rm1b))

                Rp1  = sPBC(self.path[i + 1].vdir - self.path[i].vdir)
                avgbox  = 0.5*(self.path[i + 1].get_cell() + self.path[i].get_cell())
                Rp1  = numpy.dot(Rp1,avgbox)
                dh   = self.path[i + 1].cellt - self.path[i].cellt
                Rp1b = numpy.dot(self.path[i].icell, dh)*0.5+numpy.dot(self.path[i + 1].icell, dh)*0.5
                Rp1  = numpy.vstack((Rp1,Rp1b))

                self.path[i].fsN = (vmag(Rp1) - vmag(Rm1)) * self.k * self.path[i].n
                #---------------01/26/11 to speedup by weakening spring force's convergence-----------
                #if vmag(self.path[i].fsN) < 0.01:
                    #self.path[i].fsN = 0.0
                #-------------------------------------------------------------------------------------

                # For dneb use total spring force -spring force in the grad
                # direction.
                if self.dneb:
                    self.path[i].fs = (Rp1 + Rm1) * self.k
                    self.path[i].fsperp = self.path[i].fs -                   \
                                          vproj(self.path[i].fs, self.path[i].n)
                    self.path[i].fsdneb = self.path[i].fsperp -               \
                                          vproj(self.path[i].fs, self.path[i].fPerp)

                    # New dneb where dneb force converges with (What?!)
                    if not self.dnebOrg:
                        FperpSQ = vmag(self.path[i].fPerp)
                        FsperpSQ = vmag(self.path[i].fsperp)                 
                        if FsperpSQ > 0:
                            self.path[i].fsdneb *= 2.0 / pi * atan(FperpSQ /  \
                                                                   FsperpSQ)
                
                # Not using double-nudging, so set the double-nudging spring
                # force to zero.
                else:
                    self.path[i].fsdneb = 0
                
                # The final force is the sum of these forces.    
                self.path[i].totalf = self.path[i].fsdneb + self.path[i].fsN +     \
                                 self.path[i].fPerp

                # only move the climing image
                if(self.method == 'ci' and self.onlyci): 
                    self.path[i].totalf *= 0.0
                
