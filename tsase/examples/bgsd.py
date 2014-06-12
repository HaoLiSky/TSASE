#!/usr/bin/env python

import tsase
import ase

######## load V #########################
al = tsase.calculators.al()
p = tsase.io.read_con('al.con')
p.set_calculator(al)
bmin = tsase.optimize.SDLBFGS(p)
#########################################

##### set calculator to be H landscape of the PES of p using BGSD_potential #######
p1 = tsase.io.read_con('al.con')
Hpotential = tsase.bgsd.BGSD_potential(pot = p,alpha = 5.0,beta = 0.2 ,dr=1.e-5)
p1.set_calculator(Hpotential)
###################################################################################

####### run BGSD #########
bgsd = tsase.bgsd.BGSD(p,alpha=5.0,beta=0.27,displace_atomlist=[25,26,27],displace_radius=0.0)
converged, ForceCalls = bgsd.find_saddle()

print 'Force Calls:',ForceCalls
print 'Was first order saddle point found',converged
