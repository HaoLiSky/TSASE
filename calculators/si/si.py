
import os

from ase.calculators.lammps import LAMMPS

def si(cmd = None):
    cwd = os.getcwd()
    if cmd == None:
        cmd = os.path.join(cwd, 'lmp_serial')
    os.environ['LAMMPS_COMMAND'] = cmd
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    parameters = {'pair_style': 'meam',
                  'pair_coeff': ['* * library.meam Si NULL Si']}
    files = ['library.meam']
    calc = LAMMPS(parameters = parameters, files = files)
    os.chdir(cwd)
    return calc
    
    
    
