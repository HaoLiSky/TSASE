#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Test suit for the CP2K ASE calulator with neb calculator from TSASE

http://www.cp2k.org
Author: Ole Schuett <ole.schuett@mat.ethz.ch>
"""

from __future__ import division, print_function
import os

from ase.test import NotAvailable
from ase import atoms
from ase import units
from ase.structure import molecule
from ase.calculators.cp2k import CP2K
from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from tsase import neb


"""The input for CP2K, modify accordinly"""
inp = """
 &FORCE_EVAL
   METHOD        QS
   STRESS_TENSOR  NONE
   &DFT
     BASIS_SET_FILE_NAME /opt/programs/cp2k-3.0/data/BASIS_MOLOPT
     POTENTIAL_FILE_NAME /opt/programs/cp2k-3.0/data/GTH_POTENTIALS
     UKS T
     &SCF
       MAX_ITER_LUMO     299
       MAX_SCF           200
       MAX_SCF_HISTORY   0
       MAX_DIIS          4
       LEVEL_SHIFT       0.00
       EPS_SCF           1.000E-05
       EPS_SCF_HISTORY   0.000E+00
       EPS_DIIS          1.000E-01
       SCF_GUESS         restart
       &OT                       T
         ALGORITHM               STRICT
         IRAC_DEGREE             4
         MAX_IRAC                50
         ORTHO_IRAC              CHOL
         EPS_IRAC                1.0000E-10
         EPS_IRAC_QUICK_EXIT     1.0000E-05
         EPS_IRAC_SWITCH         1.0000E-02
         ON_THE_FLY_LOC          F
         MINIMIZER               DIIS
         SAFE_DIIS               T
         N_DIIS                  5
         LINESEARCH              2PNT
         STEPSIZE                1.5000E-01
         GOLD_TARGET             1.0000E-02
         PRECONDITIONER          FULL_SINGLE_INVERSE
         PRECOND_SOLVER          DEFAULT
         ENERGY_GAP              2.0000E-01
         EPS_TAYLOR              1.0000E-16
         MAX_TAYLOR              4
         ROTATION                F
         ENERGIES                F
       &END OT
       &OUTER_SCF                T
         TYPE                    NONE
         OPTIMIZER               NONE
         BISECT_TRUST_COUNT      10
         EPS_SCF                 1.000E-05
         DIIS_BUFFER_LENGTH      3
         EXTRAPOLATION_ORDER     3
         MAX_SCF                 5
         STEP_SIZE               5.000E-01
       &END OUTER_SCF
     &END SCF
     &QS
       EPS_DEFAULT          1.000E-10
       EXTRAPOLATION        ASPC
       EXTRAPOLATION_ORDER  3
       METHOD               GPW
     &END QS
     &MGRID
       NGRIDS               5
       CUTOFF               2.800E+02
       REL_CUTOFF           3.000E+01
     &END MGRID
     &XC
       DENSITY_CUTOFF               1.0000E-10
       GRADIENT_CUTOFF              1.0000E-10
       DENSITY_SMOOTH_CUTOFF_RANGE  0.0000E+00
       TAU_CUTOFF                   1.0000E-10
       FUNCTIONAL_ROUTINE           NEW
       &XC_GRID
         XC_SMOOTH_RHO              NONE
         XC_DERIV                   PW
         USE_FINER_GRID             F
       &END XC_GRID
      &XC_FUNCTIONAL               NO_SHORTCUT
         &BECKE88                   T
           SCALE_X                  1.0000E+00
         &END BECKE88
         &LYP                       T
           SCALE_C                  1.0000E+00
         &END LYP
       &END XC_FUNCTIONAL
       &XC_POTENTIAL
         ENERGY                     NONE
       &END XC_POTENTIAL
      &vdW_POTENTIAL
         DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          REFERENCE_FUNCTIONAL BLYP
          CALCULATE_C9_TERM .TRUE.
          PARAMETER_FILE_NAME /opt/programs/cp2k-3.0/data/dftd3.dat
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END vdW_POTENTIAL
     &END XC
     &POISSON
       POISSON_SOLVER  periodic
       PERIODIC        xyz
     &END POISSON
   &END DFT

   &SUBSYS
     &KIND O 
       BASIS_SET DZVP-MOLOPT-SR-GTH-q6
       POTENTIAL GTH-BLYP-q6
     &END KIND

      &KIND H
       BASIS_SET DZVP-MOLOPT-SR-GTH-q1
       POTENTIAL GTH-BLYP-q1
     &END KIND

      &KIND N
       BASIS_SET DZVP-MOLOPT-SR-GTH-q5
       POTENTIAL GTH-BLYP-q5
     &END KIND

   &END SUBSYS
 &END FORCE_EVAL
"""

def main():
    """add the following line in .bash_profile and source the .bash_profile to
       define the varialbe "ASE_CP2K_COMMAND": 
       export ASE_CP2K_COMMAND="mpiexec.hydra -n NCORES cp2k_shell.popt"
       where NCORES is the total number of cores used for calculation.
       mpiexec.hydra is the mpi command used for the parrallel running.
    """
    if "ASE_CP2K_COMMAND" not in os.environ:
        raise NotAvailable('$ASE_CP2K_COMMAND not defined')

    # Basically, the entire CP2K input is passed in explicitly.
    # Disable ASE's input generation by setting everything to None.
    # ASE should only add the CELL and the COORD section.
    calc = CP2K(basis_set=None,
                basis_set_file=None,
                max_scf=None,
                cutoff=None,
                force_eval_method=None,
                potential_file=None,
                poisson_solver=None,
                pseudo_potential=None,
                stress_tensor=False,
                xc=None,
                label='N2O4-NH4', inp=inp)

"""Read the reactant coordinates"""
    p1 = read('reactant.xyz',index=0,format='xyz')
    p1.set_cell([[20,0,0],[0,20,0],[0,0,20]],scale_atoms=False,fix=None)
    p1.set_pbc((True, True, True))
    p1.set_calculator(calc)

"""Read the product coordinates"""
    p2 = read('product.xyz',index=0,format='xyz')
    p2.set_cell([[20,0,0],[0,20,0],[0,0,20]],scale_atoms=False,fix=None)
    p2.set_pbc((True, True, True))
    p2.set_calculator(calc)

    nim = 7  # number of images, including end points
    band = neb.ssneb(p1, p2, numImages = nim, method = 'ci', ss=False)

# to restart, uncomment the following lines which read the previous optimized images into the band
#    for i in range(1,nim-1):
#        filename = str(i)+'.CON'
#        b = read(filename,format='vasp')
#        band.path[i].set_positions(b.get_positions())
#        band.path[i].set_cell(b.get_cell())

    opt = neb.qm_ssneb(band, maxmove =0.2, dt=0.1)
#Uncomment to use fire optimization algorithm
#    opt = neb.fire_ssneb(band, maxmove =0.05, dtmax = 0.03, dt=0.03)
    opt.minimize(forceConverged=0.02, maxIterations = 200)


main()
# EOF

