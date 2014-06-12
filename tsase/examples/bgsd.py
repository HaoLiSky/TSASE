#!/usr/bin/env python

import tsase

######## load V #########################
al = tsase.calculators.al()
p = tsase.io.read_con('al.con')
p.set_calculator(al)

#### minimize p so the system is at a local minimum ##########
bmin = tsase.optimize.SDLBFGS(p)
bmin.run()

##### set calculator to be H landscape using the PES of p and the BGSD_potential #######
p1 = tsase.io.read_con('al.con')
Hpotential = tsase.bgsd.BGSD_potential(pot = p,alpha = 5.0,beta = 0.2 ,dr=1.e-5)
p1.set_calculator(Hpotential)

####### run BGSD ##########################
bgsd = tsase.bgsd.BGSD(p,alpha=5.0,beta=0.27,displace_atomlist=[25,26,27],displace_radius=0.0)

"""

The function 'find saddle' samples the beta isosurface, minimizes H and the gradient squared landscape, 
and checks to see if a 1st order saddle point connected to the reactant state is found.  
It returns two parameters. The first parameter is True if a first order connected saddle point is found and 
False otherwise.  The second parameters is the number force calls required. 
The script produces a file of called 'SP.con' which is the critical point found in the optimization.  

"""

converged, ForceCalls = bgsd.find_saddle()

###########################################
print 'Force Calls:',ForceCalls
print 'Was first order saddle point found',converged
