'''
ssneb mimizer superclass
'''

from util import vmag, sPBC
from ase import io
from numpy import dot,sqrt,vdot

class minimizer_ssneb:
    '''
    Neb minimizer superclass
    '''

    def __init__(self, band): 
        
        self.band = band

    def minimize(self, forceConverged = 0.01, maxIterations = 1000):
        '''
        Minimize the neb
            forceConverged  - stopping criterion; magnitue of the force vector
            maxForceCalls   - maximum number of force calls allowed
            maxIterations   - maximum number of iterations allowed
        '''
        fMax = 1e300
        iterations = 0
	print "Iteration       Total Force       Perp Force        MaxU       MaxI    Stress on CI    "
	print "---------------------------------------------------------------------------------------"
	feout = open('fe.out','a')
	feout.write('Iteration       Total Force       Perp Force        MaxU       MaxI    Stress on CI     \n')
	feout.write('------------------------------------------------------------  \n')
        while fMax > forceConverged and iterations < maxIterations:
            self.step()
            fMax = 0.0
            fPMax = 0.0
            for i in range(1, self.band.numImages - 1):
                fi  = vmag(self.band.path[i].totalf)
                fPi = vmag(self.band.path[i].fPerp)
	        #fi  = np.max(abs(self.band.path[i].totalf))
	        #fPi = np.max(abs(self.band.path[i].fPerp))/self.band.jacobian
                if fi > fMax:
                    fMax = fi
                if fPi > fPMax:
                    fPMax = fPi
		io.write(str(i)+'.CON',self.band.path[i],format='vasp')

            maxi=self.band.Umaxi
            fci =self.band.path[maxi].st 
            fci =vmag(fci)
            #fci =np.max(abs(fci))/self.band.jacobian
	    output = str(iterations+1)+'     '+str(fMax)+'     '+str(fPMax)+'     ' \
		     +str(self.band.Umax-self.band.path[0].u)+'     '+str(self.band.Umaxi) \
		     +'     '+str(fci) + '    '+str(self.dt)
	    print "-------------------------SSNEB------------------------------"
            print output
	    feout.write(output+'\n')

            iterations += 1

        # write data for neb.dat
	print "-----------------------SSNEB Finished------------------------------"
	print "Image    ReCoords      E      RealForce      Image"
	feout.write("-----------------------SSNEB Finished------------------------------\n")
	feout.write("Image    ReCoords      E      RealForce      Image \n")
	for i in range(self.band.numImages):
	    if i==0:
		 Rm1 = 0.0
		 R20 = 0.0
		 realtotalf = 0.0
	    else:
		 Rm1  = sPBC(self.band.path[i - 1].vdir - self.band.path[i].vdir)
		 avgb = 0.5*(self.band.path[i - 1].get_cell() + self.band.path[i].get_cell())
		 Rm1  = dot(Rm1,avgb) 
		 dh   = self.band.path[i - 1].cellt - self.band.path[i].cellt
		 Rm1b = dot(self.band.path[i].icell, dh)
		 Rm1  = sqrt(vdot(Rm1,Rm1)+vdot(Rm1b,Rm1b))
		 if i==self.band.numImages-1:
		      realtotalf = 0
		 else:
		      realtotalf = vdot(self.band.path[i].realtf,self.band.path[i].n)
	    R20 += Rm1
	    print "%3i %13.6f %13.6f %13.6f %3i" % (i,float(R20),float(self.band.path[i].u-self.band.path[0].u),float(realtotalf),i)
	    feout.write( "%3i %13.6f %13.6f %13.6f %3i %s" % (i,float(R20),float(self.band.path[i].u-self.band.path[0].u),float(realtotalf),i,'\n'))

	feout.close()
	    
	    
		
    








