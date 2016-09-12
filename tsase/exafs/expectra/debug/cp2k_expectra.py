#!/usr/bin/env python
import sys
sys.path.append("..")
import os

from ase.test import NotAvailable
from ase.io import read
from ase.optimize.lbfgs import LBFGS
from expectra.aselite import Atoms
from expectra.cal_exafs import Expectra
from expectra.io import read_xdatcar, read_con, read_chi
#from ase import *
from expectra.basin import BasinHopping
from ase.calculators.cp2k import CP2K

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
       EPS_SCF           1.000E-04
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
         EPS_SCF                 1.000E-04
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
#      &vdW_POTENTIAL
#         DISPERSION_FUNCTIONAL PAIR_POTENTIAL
#        &PAIR_POTENTIAL
#          TYPE DFTD3
#          REFERENCE_FUNCTIONAL BLYP
#          CALCULATE_C9_TERM .TRUE.
#          PARAMETER_FILE_NAME /opt/programs/cp2k-3.0/data/dftd3.dat
#          R_CUTOFF 15.0
#        &END PAIR_POTENTIAL
#      &END vdW_POTENTIAL
     &END XC
     &POISSON
       POISSON_SOLVER  periodic
       PERIODIC        xyz
     &END POISSON
   &END DFT

   &SUBSYS
     &KIND Au
       BASIS_SET SZV-MOLOPT-SR-GTH-q11
       POTENTIAL GTH-PBE-q11
     &END KIND

   &END SUBSYS
 &END FORCE_EVAL
"""

def main():
    if "ASE_CP2K_COMMAND" not in os.environ:
        raise NotAvailable('$ASE_CP2K_COMMAND not defined')

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
                label='Au55', inp=inp)

    p1 = read('geometry.xyz',index=0,format='xyz')
    p1.set_cell([[20,0,0],[0,20,0],[0,0,20]],scale_atoms=False,fix=None)
    p1.set_pbc((True, True, True))
#    filename='CONTCAR'
#    p1 = read_xdatcar(filename)
#    print(p1[0].get_positions())
#    print "p1 type in test.py"
#    print(type(p1))
    exafs_calc = Expectra(kmax = 10.0)
    bh = BasinHopping(atoms=p1,
                      opt_calculator = calc,
                      exafs_calculator=exafs_calc,
#                      temperature=300 * kB,
                      dr=0.5,
                      optimizer=LBFGS,
                      fmax=0.1)
    bh.run(2)
#    p1[0].set_calculator(exafs_calc)
#    chi_devia = p1[0].get_potential_energy()
#    print(chi_devia)
if __name__ == '__main__':
    main()
    
