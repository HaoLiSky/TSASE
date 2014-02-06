import numpy
numpy.seterr(all='raise')
import calculators
import io
import data
import neb
import kdb
import structure
import constraints
import md
import svm
import optimize
#import mc

def interpolate(trajectory, intermediate_frames = 8):
    interpolated = []
    for i in range(len(trajectory)-1):
        interpolated.append(trajectory[i].copy())
        for j in range(0, intermediate_frames):
            temp = trajectory[i].copy()
            temp.positions = trajectory[i].positions + (trajectory[i+1].positions-trajectory[i].positions) * ((j+1)/float(intermediate_frames+1))
            interpolated.append(temp)
    interpolated.append(trajectory[-1])
    return interpolated
    
    
