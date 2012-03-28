
import os

from ase.calculators.lammps import LAMMPS

def lisi(cmd = 'lmp_serial'):
    os.environ['LAMMPS_COMMAND'] = cmd
    cwd = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    parameters = {'pair_style': 'meam',
                  'pair_coeff': ['* * %s Li Si %s Li Si' % ('library.meam', 'LiSi.meam')]}
    files = ['LiSi.meam', 'library.meam']
    calc = LAMMPS(parameters = parameters, files = files)
    os.chdir(cwd)
    return calc
    
    
    
