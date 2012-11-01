'''
The lepspho (LEPS plus harmonic oscillator) potential module plus two Gaussians.
'''

import numpy 
import lepspho_

class lepspho:
    '''
    The lepspho (LEPS plus harmonic oscillator) potential class.
    '''
        
    def __init__(self):
        '''
        Initialize the lepspho (LEPS plus harmonic oscillator) potential 
        module.
        '''
        self.atoms = None
        self.u = None
        self.f = None

    def calculate(self):
        '''
        Perform a force call with the lepspho (LEPS plus harmonic oscillator) 
        potential.
        Parameters:
            p:  tsse.point object.
        '''
        ra = numpy.zeros(6,'d')
        fa = numpy.zeros(6,'d')
        uRet = numpy.array([0],'d')
        rtmp  = self.atoms.get_positions()
        ra[0] = rtmp[0][0]
        ra[3] = rtmp[0][1]
        lepspho_.force(ra,fa,uRet)
        self.f = rtmp * 0.0
        self.f[0][0] = fa[0]
        self.f[0][1] = fa[3]
        self.u = uRet[0]
        ugauss, fgauss = self.gaussianpart(ra)
        self.u += ugauss
        self.f += fgauss

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

    def gaussianpart(self, ra):      
        A1  = -1.5
        x01 = 2.02083
        y01 = -0.172881
        sigmax1 = 0.1 
        sigmay1 = 0.35
        A2  = 6.0
        x02 = 0.8
        y02 = 2.0
        sigmax2 = 0.5
        sigmay2 = 0.7
        rx  = ra[0]
        ry  = ra[3]
        ugauss1 = self.gaussian_u(rx, ry, A1, x01, y01, sigmax1, sigmay1)
        ugauss2 = self.gaussian_u(rx, ry, A2, x02, y02, sigmax2, sigmay2)
        ugauss  = ugauss1 + ugauss2
        fgauss1 = -self.gaussian_f(rx, ry, A1, x01, y01, sigmax1, sigmay1)
        fgauss2 = -self.gaussian_f(rx, ry, A2, x02, y02, sigmax2, sigmay2)
        fgauss  = fgauss1 + fgauss2
        return ugauss, fgauss
        
    def gaussian_u(self, rx, ry, A, x0, y0, sigmax, sigmay):
        gxy = A * numpy.exp(-(rx - x0)**2 / (2 * sigmax**2)) \
                * numpy.exp(-(ry - y0)**2 / (2 * sigmay**2))
        return gxy
 
    def gaussian_f(self, rx, ry, A, x0, y0, sigmax, sigmay):
        gfx = A * numpy.exp(-(rx - x0)**2 / (2 * sigmax**2)) \
                * numpy.exp(-(ry - y0)**2 / (2 * sigmay**2)) \
                * (-rx + x0) / sigmax**2
        gfy = A * numpy.exp(-(rx - x0)**2 / (2 * sigmax**2)) \
                * numpy.exp(-(ry - y0)**2 / (2 * sigmay**2)) \
                * (-ry + y0) / sigmay**2
        return numpy.array([[gfx, gfy, 0]])
