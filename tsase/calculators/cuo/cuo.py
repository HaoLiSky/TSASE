
import os

from tsase.calculators.lammps_ext import LAMMPS
#from ase.calculators.lammps import LAMMPS

def cuo(cmd = None, generation=3, tmp_dir = None, tolqeq = 0.01):
    cwd = os.getcwd()
    if cmd == None:
        cmd = os.path.join(cwd, 'lmp_serial')
    os.environ['LAMMPS_COMMAND'] = cmd
    lmptsase = os.path.dirname(os.path.abspath(__file__))
    os.chdir(lmptsase)
    if generation==2:
        pstyle = 'comb'
        files = ['ffield.comb']
    elif generation == 3:
        pstyle = 'comb3 polar_off'
        files = ['ffield.comb3']
    parameters = {'pair_style': pstyle,
                  'pair_coeff': ['* * %s Cu O' % (files[0])],
                  'atom_style':'charge',
                  'mass':['1 64','2 16'],
                  #'fix': '1 all qeq/comb 1 ' + str(tolqeq) + ' file log_qeq'}
                  'fix': '1 all qeq/comb 1.0 ' + str(tolqeq) + ' nevery 1 verlet maxloop 200'}
    if tmp_dir == None:
        calc = LAMMPS(parameters = parameters, files = files)
    else:
        tmp_dir = cwd+'/'+tmp_dir
        calc = LAMMPS(parameters = parameters, files = files, tmp_dir = tmp_dir)
    os.chdir(cwd)
    return calc
    
    
    
