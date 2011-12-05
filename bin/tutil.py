
import sys
import pickle

import ase
import tsase

def pipein():
    return pickle.loads(''.join(sys.stdin.readlines()))

def pipeout(a):
    print pickle.dumps(a, 0)

def getpot(potstr):
    try:
        exec("pot = ase.calculators.%s()" % potstr)
    except:
        exec("pot = tsase.calculators.%s()" % potstr)
    return pot
