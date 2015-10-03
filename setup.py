#!/usr/bin/env python
from distutils.core import setup
from distutils.command.install import install as DistutilsInstall
import subprocess
from sys import exit

class CustomInstall(DistutilsInstall):
        def run(self):
            subprocess.check_call('cd tsase; make', shell=True)
            DistutilsInstall.run(self)

setup(name='tsase',
      version='1.0',
      description='Library based upon ASE for transition state theory calculations.',
      author='Henkelman Research Group',
      author_email='henkelman@utexas.edu',
      url='http://www.henkelmanlab.org',
      packages=['tsase'],
      cmdclass={'install' : CustomInstall},
     )
