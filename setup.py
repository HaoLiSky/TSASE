#!/usr/bin/env python
from distutils.core import setup
from distutils.command.install import install as DistutilsInstall
import subprocess
from sys import exit

import os

class CustomInstall(DistutilsInstall):
        def run(self):
            subprocess.check_call('cd tsase; make', shell=True)
            DistutilsInstall.run(self)

packages = []
for dirname, dirnames, filenames in os.walk('tsase'):
        if '__init__.py' in filenames:
            packages.append(dirname.replace('/', '.'))

package_dir = {'tsase': 'tsase'}

scripts = ['bin/xyz']

package_data = {'tsase': ['xyz/xyz.glade',
                          'xyz/xyz.help',
                          'xyz/*.png',
                          'calculators/al/al_.so',
                          'calculators/cuo/ffield.comb',
                          'calculators/cuo/ffield.comb3',
                          'calculators/cuo/in.lammps',
                          'calculators/lepspho/lepspho_.so',
                          'calculators/lisi/LiSi.meam',
                          'calculators/lisi/in.lammps',
                          'calculators/lisi/library.meam',
                          'calculators/lj/lj_.so',
                          'calculators/mo/Mo.set',
                          'calculators/mo/in.lammps',
                          'calculators/morse/morse_.so',
                          'calculators/si/Si.meam',
                          'calculators/si/library.meam',
                          'calculators/w/W.set',
                          'calculators/w/in.lammps']}

setup(name='tsase',
      version='1.0',
      description='Library based upon ASE for transition state theory calculations.',
      author='Henkelman Research Group',
      author_email='henkelman@utexas.edu',
      url='http://www.henkelmanlab.org',
#      packages=['tsase'],
      scripts=scripts,
      packages=packages,
      package_data=package_data,
      cmdclass={'install' : CustomInstall},
     )
