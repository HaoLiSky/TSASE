import mpi4py.MPI

import argparse
import sys

import numpy

from expectra.exafs import exafs_first_shell, exafs_multiple_scattering
from expectra.io import read_xdatcar, read_con, read_chi
from expectra.feff import load_chi_dat

COMM_WORLD = mpi4py.MPI.COMM_WORLD

def mpiexcepthook(type, value, traceback):
    sys.__excepthook__(type, value, traceback)
    sys.stderr.write("exception occured on rank %i\n" % COMM_WORLD.rank)
    COMM_WORLD.Abort()
sys.excepthook = mpiexcepthook

class expectra:
    def __int__(self, kmin, kmax, trajectory, neighbor_cutoff = 6.0, s02 =
                0.89, energy_shift = 3.4, edge = 'L3', absorber = 'Au', skip =
                0, every = 1,multiple_scattering = True,
                exp_chi_file='chi_exp.dat', chi_deviation = '100')
    """
    The expectra constructor:
        kmin, kmax ...........define the k window you are interested on. It is
        suggested to set them as the values appeared in experimental data
        trajectory......coordinates or trajectories
        chi_deviaiton.........used to store the calcuation deviation between
        experimental and calculated EXAFS spectra
    """
    self.kmin = kmin
    self.kmax = kmax
    self.trajectory = trajectory
    self.neighbor_cutoff = neighbor_cutoff
    self.s02 = s02
    self.energy_shift =energy_shift
    self.edge = edge
    self.absorber = absorber
    self.skip = skip
    self.every = every

#load experimental chi data
    try:
        k_exp, chi_exp = read_chi(exp_chi_file) 
    except:
        k_exp, chi_exp = load_chi_dat(exp_chi_file)

#    last_index = len(k_exp) - 1
#    kmin = k_exp[0]
#    kmax = k_exp[last_index]
    dk   = k_exp[1] - k_exp[0]


#    if args.ignore_elements:
#        args.ignore_elements = args.ignore_elements.split(',')

    #calculate EXAFS

    trajectory = COMM_WORLD.bcast(self.trajectory)

    self.absorber = get_default_absorber(trajectory[0], self)

    k, chi = exafs_trajectory(self, trajectory)
    
    #interpolate chi values based on k values provided in experimental data
    k, chi = rescale_chi_calc(k_exp, k, chi)

    #cut a range of (k,chi) used for exp_calc_deviation cacluation 
    window = hanning_window_origin(k_exp, kmin, kmax, dk)
    chi_exp *= window
    window = hanning_window_origin(k_exp, kmin, kmax, dk)
    chi *= window

    save_result(k, chi)

    self.chi_diviation = exp_calc_deviation(chi_exp, chi)

    #need to modify it to save a history of (chi, k)
    def save_result(k, chi):
        if COMM_WORLD.rank != 0: return
        print 'saving result to chi.dat'
        f = open('chi.dat', 'w')
        for i in xrange(len(k)):
            f.write("%6.3f %16.8e\n" % (k[i], chi[i]))
        f.close()
    
    def get_default_absorber(atoms, args):
        symbols = set(atoms.get_chemical_symbols())
        if args.absorber:
            if args.absorber not in symbols:
                print 'ERROR: --absorber %s is not in the system' % args.absorber
                sys.exit(2)
            else:
                return args.absorber
        if args.ignore_elements:
            symbols -= set(args.ignore_elements)
        if len(symbols) == 1:
            return list(symbols)[0]
        else:
            print 'ERROR: must specify --absorber if more than one chemical specie'
            sys.exit(2)
 
    def exafs_trajectory(args, trajectory):
        if args.multiple_scattering:
            k, chi = exafs_multiple_scattering(args.S02, args.energy_shift, 
                    args.absorber, args.ignore_elements, args.edge, args.rmax, 
                    trajectory)
        elif args.first_shell:
            k, chi = exafs_first_shell(args.S02, args.energy_shift, 
                    args.absorber, args.ignore_elements, args.edge, 
                    args.neighbor_cutoff, trajectory)
    
        return k, chi

    #linearly interpolate chi values based on experimental k value
    def rescale_chi_calc(k_exp, chi_calc, k_cacl)
  
        #store original chi data
        k_temp   = k_cacl
        chi_temp = chi_calc
    
        #reset chi_calc based on k_exp
        #tell if k_exp starts from a smaller value
        try:
            result = compareValue(k_exp[0],k_cacl[0])
        except MyValidationError as exception:
            print exception.message
        
        for i in range(0,len(k_exp))
            for j in range(0,len(k_cacl))
                if k_cacl[j] < k_exp[i] and k_exp[i] < k_cacl[j+1]:
                    chi_calc[i] = numpy.interp(k_exp[i],
                                               [k_calc[j],k_cacl[j+1]],
                                               [chi_calc[j],chi_calc[j+1]])
                elif k_exp[i] == k_cacl[j]:
                    chi_calc[i] = chi_temp[j]
    
        k_cacl = k_exp
        return k_cacl, chi_calc
        
#calculate the deviation of theoretical EXAFS from experimental EXAFS
    def exp_calc_deviation(chi_exp,chi_calc)
        chi_devi = 0
        for i in range(0, len(chi_exp)):
            chi_devi = chi_devi + (chi_exp - chi_theory)**2

            if chi_exp > 0.00:
                numb += 1
        return chi_devi/numb
   
"""calculate the virtual potential which is defined as
U=alpha*chi_devi+(1-alpha)*System_energy. Alpha is the weight"""

#    def virtualPot(chi_devi,cluster_E,alpha):
#        virtual_U = alpha*chi_devi + *cluster_E
#        return virtual_U

