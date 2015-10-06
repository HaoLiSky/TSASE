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

setup(name='tsase',
      version='1.0',
      description='Library based upon ASE for transition state theory calculations.',
      author='Henkelman Research Group',
      author_email='henkelman@utexas.edu',
      url='http://www.henkelmanlab.org',
#      packages=['tsase'],
      scripts=scripts,
      packages=packages,
      cmdclass={'install' : CustomInstall},
     )
