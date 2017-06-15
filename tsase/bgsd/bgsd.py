
import tsase
import ase
import numpy
import random

### Define a potential that converts another potential into the H landscape ############
########################################################################################

class BGSD_potential:

    def __init__(self,pot,alpha=1.0,beta=1.0,dr=1e-5):
        
        """
        Parameters:

        pot - an atoms object from ASE which is the local minima in which you want to find all saddle points and calculator
        alpha - spring constant for bias term in BGSD
        beta - isosurface energy of bias term in BGSD
        dr - finite difference size for gradient of H
        
        """

        self.alpha = alpha
        self.beta = beta
        self.pot = pot
        self.dr = dr
        self.atoms = None
        self.u = None
        self.f = None
        self.g2 = None

    def calculate(self):
        self.u, self.g2, self.f = self.potforce()

    def get_gradient_squared(self, atoms=None, force_consistent=False):
        self.calculate()
        return self.g2

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
        if self.f == None or self.u == None or atoms == None or self.g2 == None:
            return True
        return False

    def set_atoms(self, atoms):
        pass

    def potforce(self):
        u,g2  = self.u1d()
        f  = self.f1d()
        return u, g2, f

    def u1d(self):
        self.pot.set_positions(self.atoms.get_positions())
        f = self.pot.get_forces()
        g2 = numpy.vdot(f,f)
        u = self.pot.get_potential_energy()
        har = self.alpha * (u - self.beta) * (u - self.beta)
        toteng = g2 + har
        return toteng,g2

    def f1d(self):
        self.pot.set_positions(self.atoms.get_positions())
        f = self.pot.get_forces()
        magf =numpy.sqrt(numpy.vdot(f,f))
        fnorm = f/magf
        u = self.pot.get_potential_energy()
        pos = self.pot.get_positions()
        gradhar = 2*self.alpha*(u - self.beta)*f
        self.pot.set_positions(pos - self.dr*fnorm)
        gradf = 2*(self.pot.get_forces() - f)*magf/self.dr
        self.pot.set_positions(pos)
        totf = gradf +gradhar
        return totf

##### Samples isosurface, min H and V2, and checks for codim1 connected saddle  ############

class BGSD:

    def __init__(self,p,alpha=5.0,beta=0.0,dr=1.0e-5,numstep=1000,stepsize=0.05,kT_MC=0.01,k=5.0,displace_atomlist=[],displace_radius=3.3,displace_all_listed=True,CC1=0.01,CC2=0.001*0.001,CC3 = 0.0001,maxstepsize=0.2,maxnumstep=1000,memory=100,eigcheck = 'dimer'):

        """
        Parameters:
        
        BGSD: p,alpha,beta, and dr analagous to parameters defined in BGSD_potential
        
        Initial configuration/Monte Carlo sampling of isosurface:
        
        numstep - number of Monte carlo steps
        stepsize - monte carlo step size
        kT_MC - temperature of MC sampling
        k - spring constant for isosurface sampling 

        Biasing initial configuration:
        
        displacement_atomlist - index number of atoms displaced in monte carlo sampling
        displace_radius - displace atoms within this distance from atomlist atoms
        displace_all_listed - if True then uses all atoms is atomlist else it randomly chooses one atom

        SDLBFGS optimizer options:

        CC1 - convergence criteria for H; L2 norm of force
        CC2 - energy convergence crtiria of gradient squared landscape
        CC3 - L2 norm criteria for gradient squared landscape (additional CC for IP) 
        maxstepsize - max step for SDLBFGS optimizer
        maxnumstep - maximum number of steps in optimization
        memory - number of terms in memory for LBFGS

        Options for minmode and second lowest eigenvalue checks:

        eigcheck - option to use hessian or dimer 

        """
    
        self.p = p
        self.pot = self.p.copy()
        self.pot.set_calculator(self.p.get_calculator())
        self.alpha = alpha
        self.beta = beta + self.p.get_potential_energy()
        self.stepsize = stepsize
        self.kT_MC = kT_MC
        self.k = k
        self.displace_atomlist = displace_atomlist
        self.displace_radius = displace_radius
        self.displace_all_listed = displace_all_listed
        self.numstep = numstep
        self.CC1 = CC1
        self.CC2 = CC2
        self.CC3 = CC3
        self.maxstepsize = maxstepsize
        self.maxnumstep = maxnumstep
        self.min_pos = self.p.get_positions()
        self.memory = memory
        self.dr = dr
        self.eigcheck = eigcheck

###### calculates hessian ######################################
    def hessian(self,numfree,clist):
        hessian = numpy.zeros((numfree*3,numfree*3))
        a=0
        dr=self.dr
        coord = self.p.get_positions()
        for i in range(len(coord)):
            if clist[i] == 1:
                for j in range(3):
                    f1 = []
                    f2 = []
                    coord[i][j] += dr
                    self.p.set_positions(coord)
                    f = self.p.get_forces()
                    for k in range(len(coord)):
                        if clist[k] == 1:
                            for m in range(3):
                                f1 = numpy.append(f1,f[k][m])

                    coord[i][j] -= (2*dr)
                    self.p.set_positions(coord)
                    f = self.p.get_forces()
                    for k in range(len(coord)):
                        if clist[k] == 1:
                            for m in range(3):
                                f2 = numpy.append(f2,f[k][m])
                    h = (f2-f1)/(2.0*dr)   
                    hessian[a]=h
                    coord[i][j] += dr
                    a += 1
        return hessian

####### Monte Carlo step for Isosurface ########
    def MC_step(self,atomlist):
        orig_pos = self.p.get_positions()
        u = self.p.get_potential_energy()
        har = self.k * (u - self.beta) * (u - self.beta)
        orig_eng = har
        new_pos = self.p.get_positions()
        for i in range(len(atomlist)):
            for j in range(3):
                new_pos[atomlist[i]][j] += random.uniform(-self.stepsize,self.stepsize)  #,seed=numpy.random.random_integers(10,50000))
        self.pot.set_positions( new_pos)
        u = self.pot.get_potential_energy()
        har = self.k * (u - self.beta) * (u - self.beta) 
        new_eng = har
        if new_eng <= orig_eng:
            self.p.set_positions(self.pot.get_positions())
        else:
            prob = numpy.exp(-(new_eng-orig_eng)/self.kT_MC)
            if numpy.random.random() <= prob:
                self.p.set_positions(self.pot.get_positions())

######## Get atomlist for MC simualtion  ########
    def make_atom_list_mc(self):
        if self.displace_all_listed == True:
            atomlist = self.displace_atomlist
            newlist = self.displace_atomlist
            for i in range(len(atomlist)):
                for j in range(len(self.p.get_positions())):
                    diff = self.p.get_positions()[atomlist[i]] - self.p.get_positions()[j]
                    magdiff = numpy.sqrt(numpy.vdot(diff,diff))
                    if magdiff < self.displace_radius:
                        flag = 0
                        for w in range(len(newlist)):
                            if newlist[w] == j:
                                flag = 1
                        if flag == 0:
                            newlist = numpy.append(newlist,j)
        else:
            ranint = numpy.random.random_integers(0,len(self.displace_atomlist)-1)
            atomlist = [ranint]
            newlist = [self.displace_atomlist[ranint]]
            for i in range(len(self.p.get_positions())):
                if i != newlist[0]:
                    diff = self.p.get_positions()[atomlist[0]] - self.p.get_positions()[i]
                    magdiff = numpy.sqrt(numpy.vdot(diff,diff))
                    if magdiff < self.displace_radius:
                        newlist = numpy.append(newlist,i)
        return newlist

###### monte carlo paramters / run initial MC of isosurface ##############
    def sample_iso(self):
        atomlist = self.make_atom_list_mc()
        print 'atomlist',atomlist   
        for i in range(self.numstep):
            self.MC_step(atomlist)
            
###### set up potential for BGSD with alpha and beta #######
    def min_H(self):
        bg = self.p.copy() 
        bg.set_positions(self.p.get_positions())
        bg.set_calculator(BGSD_potential(pot = self.p,alpha = self.alpha,beta = self.beta,dr=self.dr))
        bg.set_velocities(0.0*self.p.get_positions())
        ### minimize H #################################
        bmin = tsase.optimize.SDLBFGS(bg,logfile=None,maxstep=self.maxstepsize,memory=self.memory) #,trajectory='minH.traj')
        bmin.run(fmax=self.CC1,steps=self.maxnumstep,optimizer='L2')
        FC = 2*bmin.get_number_of_steps()
        self.p.set_positions(bg.get_positions())
        return FC

### minimize gradient squared ##################
    def min_V2(self):
        bg = self.p.copy() 
        bg.set_positions(self.p.get_positions())
        bg.set_calculator(BGSD_potential(pot = self.p,alpha = 0.0,beta = 1.,dr=self.dr))
        force = numpy.sqrt(numpy.vdot(bg.get_forces(),bg.get_forces()))
        bmin = tsase.optimize.SDLBFGS(bg,logfile=None,maxstep=self.maxstepsize,memory=self.memory) 
        bmin.run(fmax=self.CC3,emax=self.CC2,steps=self.maxnumstep,optimizer='bgsd')
        FC2 = 2*bmin.get_number_of_steps()
        force = numpy.sqrt(numpy.vdot(bg.get_forces(),bg.get_forces()))
        self.p.set_positions(bg.get_positions())
        return FC2,bg.get_potential_energy(),force

### check to see if gradV^2 has converged ###########
    def Convergence_check(self,eng,force):
        print 'convergence check',eng
        if eng > self.CC2:
            return False
        else:
            return True

##### check to see if we are at the reactant state ######
    def Check_reactant(self,distcutoff=0.5,min_pos = 0.0):
        dist = numpy.sqrt(numpy.vdot((min_pos - self.p.get_positions()),(min_pos - self.p.get_positions())))
        if dist < distcutoff:
            return True
        else:
            return False

##### get number of free atoms and list #####
    def free_atoms(self):
        frozenatoms =  self.p.constraints[0].index
        clist = numpy.ones(len(self.p.get_positions()))
        for i in range(len(frozenatoms)):
            clist[frozenatoms[i]] = 0
        numfree = len(self.p.get_positions()) - len(frozenatoms)
        return numfree,clist

### check number of negative eigenvalues by calculting hessian or using dimer ##############
    def eig_check(self):
        if self.eigcheck == 'hessian':
            numfree,clist = self.free_atoms()
            hes = self.hessian(numfree,clist)
            val,vec=numpy.linalg.eig(hes)
            print val
            counter = 0
            for i in range(len(val)):
                if val[i] < 0:
                    counter += 1
                    index = i
            if counter == 1:
                return vec[:,index].reshape((numfree,3)),counter
            else:
                return numpy.zeros(numpy.shape(vec[0])),counter
        elif self.eigcheck == 'dimer':
            from tsase.dimer import ssdimer
            numfree,clist = self.free_atoms()
            dimer1 = self.p.copy()
            dimer1.set_calculator(self.p.get_calculator())
            d = ssdimer.SSDimer_atoms(dimer1, rotationMax = 100, maxStep = 0.10, phi_tol= 0.01, ss = False, dT = 0.1)
            d.minmodesearch()
            if d.curvature > 0.0:
                return numpy.zeros((numfree,3)),0
            else:
                vec = d.get_mode()
                for i in range(3):
                    vec = numpy.delete(vec,len(vec)-1,axis=0)
                for i in range(len(vec)):
                    if clist[i] == 0:
                        vec[i] = [0.,0.,0.]
                vec /= numpy.sqrt(numpy.vdot(vec,vec))
                dimer = self.p.copy()
                dimer2 = self.p.copy()
                dimer2.set_calculator(self.p.get_calculator())
                dimer.set_calculator(tsase.calculators.hyperplane_potential(dimer2,negeig=vec))
                d1 = ssdimer.SSDimer_atoms(dimer, rotationMax = 300, maxStep = 0.10, phi_tol= 0.01, ss = False, dT = 0.1)
                d1.minmodesearch()
                if d1.curvature < -1e-3:
                    print 'FOUND HIGHER ORDER SADDLE WITH EIGENVALUE',d1.curvature
                    return numpy.zeros((numfree,3)),2
                else:
                    vec1 = numpy.zeros((numfree,3))
                    count = 0
                    for i in range(len(vec)):
                        if clist[i] == 1:
                            vec1[count] = vec[i]
                            count += 1
                    return vec1,1



##### check if saddle is connected to reactant state #####
    def check_connected(self,min_eng,min_pos,eigvec):
        numfree,clist = self.free_atoms()
        self.pot.set_positions(self.p.get_positions())
        direction = 0.5 * eigvec
        r = self.p.get_positions()
        count = 0
        for i in range(len(clist)):
            if clist[i] == 1:
                r[i] += direction[count]
                count += 1

        self.p.set_positions(r)
        bmin = tsase.optimize.SDLBFGS(self.p,logfile=None,maxstep=0.1)  #,trajectory='1.traj')
        bmin.run(fmax=0.01)
        FC = bmin.get_number_of_steps()
        filename = 'min.con'
        tsase.io.write_con(filename,self.p,w='w')

        dist1 = min_pos - self.p.get_positions()
        dist1 = numpy.sqrt(numpy.vdot(dist1,dist1))
        x = dist1 % self.p.get_cell()[0][0]
        if x > self.p.get_cell()[0][0]/2:
            x = self.p.get_cell()[0][0] - x
        else:
            dist1 = x
        count = 0
        for i in range(len(clist)):
            if clist[i] == 1:
                r[i] = self.pot.get_positions()[i] - direction[count]
                count += 1

        self.p.set_positions(r)
        bmin = tsase.optimize.SDLBFGS(self.p,logfile=None,maxstep=0.1) #,trajectory='2.traj')
        bmin.run(fmax=0.01)
        FC += bmin.get_number_of_steps()
        dist = min_pos - self.p.get_positions()
        filename = 'product.con'
        tsase.io.write_con(filename,self.p,w='w')

        dist = numpy.sqrt(numpy.vdot(dist,dist))
        x = dist % self.p.get_cell()[0][0]
        if x > self.p.get_cell()[0][0]/2:
            dist = self.p.get_cell()[0][0] - x
        else:
            dist = x
        self.p.set_positions(self.pot.get_positions())
        if dist1 < 0.5 and dist < 0.5:
            print 'not connected'
            return False,FC
        elif dist1 > 0.5 and dist > 0.5:
            print 'not connected'
            return False,FC
        else:
            print 'connected'
            return True,FC

########## gets initial configuration with sampling isosurface; minimize H and grad^2; check if connected and for codim1 SP

    def find_saddle(self):  ### assumes initial configuration is at the minimum 
        FC = 0
        min_eng = self.p.get_potential_energy()
        min_pos = self.p.get_positions()
        atomlist = self.make_atom_list_mc()
        self.sample_iso()
        FC = self.min_H()
        FCv2,eng,forces = self.min_V2()
        FC += FCv2
        filename = 'SP.con'
        tsase.io.write_con(filename,self.p, w='w')
        if self.Convergence_check(eng,forces) == False:
            print 'did not converge'
            return False,FC
        elif self.Check_reactant(min_pos = min_pos) == True:
            print 'converged at reactant state'
            return False,FC
        else:
            eigvec,numnegmodes = self.eig_check()
            if numnegmodes == 1:
                check,FCconnect = self.check_connected(min_eng,min_pos,eigvec)
                print 'FC connect', FCconnect
                FC += FCconnect
                if check == True:
                    print 'found co-dimension 1 saddles with barrier', self.p.get_potential_energy() - min_eng
                    return True,FC
                else:
                    print 'not connected to reactant state with barrier',self.p.get_potential_energy() - min_eng
                    return False,FC
            elif numnegmodes > 0:
                print 'found co-dimension',numnegmodes,'saddle with barrier',self.p.get_potential_energy() - min_eng
                return False,FC
            else:
                print 'found product state'
                return False,FC



