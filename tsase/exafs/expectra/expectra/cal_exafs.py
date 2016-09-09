#!/usr/bin/env python

import mpi4py.MPI
from expectra.MPI_Import import mpi_import
import numpy

with mpi_import():
    import argparse
    import sys

    from expectra.exafs import exafs_first_shell, exafs_multiple_scattering
    from expectra.io import read_xdatcar, read_con, read_chi
    from ase.calculators.calculator import Calculator, all_changes, Parameters
    from expectra.aselite import Atoms
    from expectra.feff import load_chi_dat

COMM_WORLD = mpi4py.MPI.COMM_WORLD

def mpiexcepthook(type, value, traceback):
    sys.__excepthook__(type, value, traceback)
    sys.stderr.write("exception occured on rank %i\n" % COMM_WORLD.rank)
    COMM_WORLD.Abort()
sys.excepthook = mpiexcepthook

#need to modify it to save a history of (chi, k)
def save_result(k, chi, filename):
    if COMM_WORLD.rank != 0: return
    print 'saving result to chi.dat'
    f = open(filename, 'w')
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
#    print "trajectory: "
#    print(type(trajectory))
    if args.multiple_scattering:
        k, chi = exafs_multiple_scattering(args.S02, args.energy_shift, 
                args.absorber, args.ignore_elements, args.edge, args.rmax, 
                trajectory)

    elif args.first_shell:
        k, chi = exafs_first_shell(args.S02, args.energy_shift, 
                args.absorber, args.ignore_elements, args.edge, 
                args.neighbor_cutoff, trajectory)
    
    return k, chi

#calculate the deviation of theoretical EXAFS from experimental EXAFS
def calc_deviation(chi_exp,chi_theory):
    chi_exp_array = numpy.asarray(chi_exp)
    chi_thry_array = numpy.asarray(chi_theory)
    chi_devi = numpy.sum(numpy.square(chi_exp_array - chi_thry_array))

    return chi_devi/len(chi_exp)

#linearly interpolate chi values based on k_std value
def rescale_chi_calc(k_std, chi_src, k_src, kmax):
    """
    k_std..........k values used as a standard for the rescaling
    chi_src........chi values required to be rescaled
    k_src..........k values corresponding to chi_src
    """
    k_temp = []
    chi_temp = []
    #reset chi_calc based on k_exp
    #tell if k_exp starts from a smaller value
#          try:
#          result = compareValue(k_exp[0],k_cacl[0])
#      except MyValidationError as exception:
#          print exception.message
    i = 0   
    while ( 0 <= i < len(k_std) and k_std[i] < kmax):
        for j in range(1,len(k_src)):
            if k_src[j-1] < k_std[i] and k_std[i] < k_src[j]:
                chi_temp.append(numpy.interp(k_std[i],
                                           [k_src[j-1],k_src[j]],
                                       [chi_src[j-1],chi_src[j]]))
                k_temp.append(k_std[i])

            elif k_std[i] == k_src[j-1]:
                chi_temp.append(chi_src[j-1])
                k_temp.append(k_std[i])
        i += 1
    return k_temp, chi_temp

'''
Calculator is the superclass. expectra is the subclass
'''
class Expectra(Calculator):

    implemented_properties = ['chi_deviation']

    default_parameters = dict(
        multiple_scattering = True,
        ignore_elements = None,
        neighbor_cutoff = 6.0,
        rmax = 6.0,
        S02 = 0.89,
        energy_shift = 3.4,
        edge = 'L3',
        absorber = 'Au',
        skip = 0,
        exp_chi_file = 'chi_exp.dat',
        output_file = 'chi.dat')

    def __int__(self, label='EXAFS', 
                atoms=None, kmin=0.00, kmax=10.00, chi_deviation=100, **kwargs):
        """
        The expectra constructor:
            kmin, kmax ...........define the k window you are interested on. It is
            suggested to set them as the values appeared in experimental data
            atoms.................coordinates or trajectories
            chi_deviaiton.........used to store the calcuation deviation between
            experimental and calculated EXAFS spectra
        """
        self.lable = label
        self.atoms = atoms
        self.kmin = kmin
        self.kmax = kmax
        self.chi_deviation = chi_deviation
        self.parameters = None
        self.results = None
      
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms,
                            **kwargs)
      

    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2,
        ...)."""
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
           self.reset()

    def get_potential_energy(self, atoms=None, force_consistent=False):
        print(type(atoms))
        self.calculate(atoms, 'chi_deviation')
        return self.chi_deviation

    def calculate(self, atoms=None, properties=None):

        parameters = self.parameters
        print "atoms: " 
        print(type(atoms))
#        print(atoms.get_positions())
        trajectory = []
        trajectory.append(atoms)
#        print(type(trajectory))
        trajectory = COMM_WORLD.bcast(trajectory)
#        self.absorber = get_default_absorber(trajectory, parameters)
        k, chi = exafs_trajectory(parameters, trajectory)

        #load experimental chi data
        try:
            k_exp, chi_exp = read_chi(parameters.exp_chi_file) 
        except:
            k_exp, chi_exp = load_chi_dat(parameters.exp_chi_file)

        #interpolate chi_exp values based on k values provided in calculated data
        k_exp, chi_exp = rescale_chi_calc(k, chi_exp, k_exp, parameters.kmax)
        k, chi = rescale_chi_calc(k, chi, k, parameters.kmax)

        filename1 = 'chi.dat'
        filename2 = 'rescaled_exp_chi.dat'
        save_result(k, chi, filename1)
        save_result(k_exp, chi_exp, filename2)
    
        self.chi_deviation = calc_deviation(chi_exp, chi)

