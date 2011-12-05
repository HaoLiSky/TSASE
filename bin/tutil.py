
import sys
import pickle

def pipein():
    return pickle.loads(''.join(sys.stdin.readlines()))

def pipeout(a):
    print pickle.dumps(a, 0)
