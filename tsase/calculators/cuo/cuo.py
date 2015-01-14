
import os

from tsase.calculators.lammps_ext import LAMMPS
#from ase.calculators.lammps import LAMMPS

def cuo(cmd = None, tmp_dir = None, tolqeq = 0.01):
    cwd = os.getcwd()
    if cmd == None:
        cmd = os.path.join(cwd, 'lmp_serial')
    os.environ['LAMMPS_COMMAND'] = cmd
    lmptsase = os.path.dirname(os.path.abspath(__file__))
    os.chdir(lmptsase)
    parameters = {'pair_style': 'comb',
                  'pair_coeff': ['* * %s O Cu' % ('ffield.comb')],
                  'atom_style':'charge',
                  'mass':['1 16','2 64'],
                  'fix': '1 all qeq/comb 1 ' + str(tolqeq) + ' file log_qeq'}
    files = ['ffield.comb']
    if tmp_dir == None:
        calc = LAMMPS(parameters = parameters, files = files)
    else:
        tmp_dir = cwd+'/'+tmp_dir
        calc = LAMMPS(parameters = parameters, files = files, tmp_dir = tmp_dir)
    os.chdir(cwd)
    return calc
    
    
    
