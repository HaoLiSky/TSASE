
import os

from ase.calculators.lammpsrun import LAMMPS

def w(cmd = None):
    cwd = os.getcwd()
    if cmd == None:
        cmd = os.path.join(cwd, 'lmp_serial')
    os.environ['LAMMPS_COMMAND'] = cmd
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    parameters = {'pair_style': 'eam/alloy',
                  'pair_coeff': ['* * W.set W']}
    files = ['W.set']
    calc = LAMMPS(parameters = parameters, files = files)
    os.chdir(cwd)
    return calc
    
    
    
