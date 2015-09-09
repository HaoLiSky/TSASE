from numpy import exp
import numpy 
import math

try:
    import sklearn
except:
    import scikits.learn

class svm_dynamics():
    def __init__(self,calc,selSVM,k=1.0,addgradient=True):
        
        self.atoms = None
        self.calc = calc
        self.selSVM = selSVM
        self.u = None
        self.f = None
        self.k = k
        self.addgradient = addgradient

    def set_k(self,k):
        self.k = k
        
    def get_k(self):
        return self.k

    def include_SVM_gradient(self,addgradient):
        self.addgradient = addgradient

    def gradient_included(self):
        return self.addgradient

    def decision(self,atoms):
        self.atoms = atoms.copy()
        datapoint = numpy.asarray([self.atoms.get_positions().ravel()])
        decisionf = self.selSVM.decision_function(datapoint)[0,0]
        return decisionf

    def calculate(self):
        self.calc.get_forces(self.atoms)
        self.f = self.calc.f.copy()
        if self.addgradient:
            datapoint = numpy.asarray([self.atoms.get_positions().ravel()])
            grad = self.svm_rbf_deriv(self.selSVM, datapoint) 
            gradnorm = numpy.sqrt(numpy.vdot(grad,grad))
            if (gradnorm >= 1e-10):
                grad /= gradnorm # the gradient is now normalized
                prediction = self.selSVM.decision_function(datapoint)[0,0] / gradnorm
            else:
                prediction = 0.
            self.f += -1 * self.k * prediction * numpy.reshape(grad,numpy.shape(self.atoms.get_positions()))

    def svm_rbf_deriv(self,clf,Points): # Don't forget this returns the positive of the gradient!
        import numpy as np
        try:
            from sklearn.metrics.pairwise import rbf_kernel
        except:
            from scikits.learn.metrics.pairwise import rbf_kernel
        Grad = np.zeros(Points.shape)
        noP = Points.shape[0]
        alpha = clf.dual_coef_
        d = np.ones((len(clf.support_),1))
        # loop over test points an compute the gradient individually
        for k in range(noP):
            P = Points[k,:]
            P.shape = 1,-1
            # be careful with the gamma paremeter: gamma= 1/(2sigma^2)
            if (float(sklearn.__version__) <= 0.8):
                Kmatrix = rbf_kernel(clf.support_vectors_,P,np.sqrt(1./(clf.gamma*2)))
            else: 
                Kmatrix = rbf_kernel(clf.support_vectors_,P,clf.gamma)
            kstar = alpha.T *Kmatrix    
            b = (clf.support_vectors_-np.dot(d,P))*2*clf.gamma
            Grad[k,:] = np.dot(kstar.T,b)
        Grad *= -1 # The function returns the negative of the gradient without this change
        return Grad

    def get_potential_energy(self, atoms, force_consistent=False):
        if self.calculation_required(atoms, "energy"):
            self.atoms = atoms.copy()
            self.calc.atoms = atoms.copy()
            self.calc.calculate()
        return self.calc.u

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

class ocsvm_dynamics(): 

    def __init__(self, p, psurf, clf, cdsa_cc=0.001, cdsa_maxstep=1.0, w=1.0, B=1.0, fixcm=False, dt=1.0, kT=0.1):
        self.p = p
        self.psurf = psurf    
        self.clf = clf
        self.cdsa_cc = cdsa_cc
        self.cdsa_maxstep = cdsa_maxstep
        self.w = w
        self.B = B
        self.fixcm = fixcm
        self.dt = dt
        self.kT = kT

    def svm_rbf_deriv(self,clf,Points): # Don't forget this returns the positive of the gradient!
	try:
            from sklearn.metrics.pairwise import rbf_kernel
        except:
            import scikits.learn as sklearn
            from scikits.learn.metrics.pairwise import rbf_kernel
        Grad = numpy.zeros(Points.shape)
        noP = Points.shape[0]
        alpha = clf.dual_coef_
        d = numpy.ones((len(clf.support_),1))
        # loop over test points an compute the gradient individually
        for k in range(noP):
            P = Points[k,:]
            P.shape = 1,-1
            # be careful with the gamma paremeter: gamma= 1/(2sigma^2)
            if (float(sklearn.__version__) <= 0.8):
                Kmatrix = rbf_kernel(clf.support_vectors_,P,numpy.sqrt(1./(clf.gamma*2)))
            else:
                Kmatrix = rbf_kernel(clf.support_vectors_,P,clf.gamma)
            kstar = alpha.T *Kmatrix
            b = (clf.support_vectors_-numpy.dot(d,P))*2*clf.gamma
            Grad[k,:] = numpy.dot(kstar.T,b)
        Grad *= -1 # The function returns the negative of the gradient without this change
        return Grad

    def initguess(self):
        pt = self.p.get_positions().ravel()
        svar = self.clf.support_vectors_
        max3 = numpy.zeros((1,len(svar[0])))
        for i in range(len(svar)):
            diff = pt - svar[i]
            dist = numpy.sqrt(numpy.dot(diff,diff))
            if i == 0:
                maxdist = dist
                max = svar[i]
                index = i
            else:
                if dist < maxdist:
                    maxdist = dist
                    max = svar[i]
                    index = i
        max3[0] = max
        return max3

    def NM(self):
        pt = self.p.get_positions().ravel()
        distvec = self.psurf.get_positions().ravel() - pt
        dist = numpy.sqrt(numpy.vdot(distvec,distvec))
        newtonstep = 1.0
        while newtonstep > self.cdsa_cc:
            distvec = self.psurf.get_positions().ravel() - pt
            dist = numpy.sqrt(numpy.vdot(distvec,distvec))
            ndistvec = distvec/dist
            deriv = self.svm_rbf_deriv(self.clf,numpy.asarray([self.psurf.get_positions().ravel()]))
            ss = numpy.vdot(deriv,ndistvec)
            df = self.clf.decision_function(numpy.asarray([self.psurf.get_positions().ravel()]))
            if ss != 0: #### prob should be not equal to instead of greater than
                newtonstep = df/ss
            else:
                newtonstep = 0.01
            if math.fabs(newtonstep) > self.cdsa_maxstep:
                        newtonstep = self.cdsa_maxstep
            if math.fabs(newtonstep) > self.cdsa_cc:
                ndistvec = numpy.reshape(ndistvec,numpy.shape(self.psurf.get_positions()))
                self.psurf.set_positions(self.psurf.get_positions() + newtonstep*ndistvec)
            newtonstep = math.fabs(newtonstep)
        return deriv

    def find_dist_surf(self):
        pt = self.p.get_positions().ravel()
        svar = self.initguess()
        sv = numpy.reshape(svar[0],numpy.shape(self.p.get_positions()))
        self.psurf.set_positions(sv)
        a = self.NM()
        magp = 1
        count = 0
        while magp >  self.cdsa_cc:
            anorm = numpy.sqrt(numpy.vdot(a,a))
            if anorm > 0:
                a /= anorm
            b = pt - self.psurf.get_positions().ravel()
            perp = b - numpy.vdot(a,b)*a
            magperp = numpy.sqrt(numpy.vdot(perp,perp))
            if magperp > self.cdsa_maxstep:
                perp = (1.0/magperp) * perp
                perp = numpy.reshape(perp,numpy.shape(self.psurf.get_positions()))
                self.psurf.set_positions( self.psurf.get_positions() + self.cdsa_maxstep * perp)
            else:
                perp = numpy.reshape(perp,numpy.shape(self.psurf.get_positions()))
                self.psurf.set_positions( self.psurf.get_positions() + perp)
            a = self.NM()
            magp = numpy.sqrt(numpy.vdot(perp,perp))
            cdist = numpy.sqrt(numpy.vdot(self.psurf.get_positions().ravel()-pt,self.psurf.get_positions().ravel()-pt))
            magp = magp/cdist
            count += 1
        diff = self.psurf.get_positions().ravel() - pt
        mdiff = self.psurf.get_positions().ravel() - pt
        mdist = numpy.sqrt(numpy.vdot(diff,diff))
        return mdist,mdiff

    def bias_force(self):
        maxvalue = self.w*self.w*self.w/6
        Bn = self.B/maxvalue
        dist,distvec = self.find_dist_surf()
     #   print dist
        if dist < self.w:
            c = Bn*(self.w-dist)
            umforce = c*distvec
            bias = Bn*(self.w*dist*dist/2 - dist*dist*dist/3)
            return umforce,bias  
        else:
            umforce = 0.0*self.p.get_positions()
            return umforce,self.B  

    def MD_step_bias(self,umforce):
        # Note this is only the Verlet algorithm; we need to apply a thermostat separately
        if self.clf.decision_function(numpy.asarray([self.p.get_positions().ravel()])) > 0:
            f = self.p.get_forces() + umforce
        else:
            f = self.p.get_forces()
        m = self.p.get_momenta()
        m += 0.5 * self.dt * f
        if self.fixcm == True: 
                psum = m.sum(axis=0) / float(len(self.p))
                m = m - psum
        self.p.set_positions(self.p.get_positions() + self.dt * m / self.p.get_masses()[:,numpy.newaxis])

        self.p.set_momenta(m)
        if self.clf.decision_function(numpy.asarray([self.p.get_positions().ravel()])) > 0:
            umforce,bias = self.bias_force()
            umforce = numpy.reshape(umforce,numpy.shape(self.p.get_forces()))
            f = self.p.get_forces() + umforce
        else:
            f = self.p.get_forces()
            umforce,bias = 0.0*self.p.get_positions(),0.0
        self.p.set_momenta(self.p.get_momenta() + 0.5 * self.dt * f)
        boost = numpy.exp(bias/self.kT)
        time = self.dt*boost
        return time,boost,umforce,bias

