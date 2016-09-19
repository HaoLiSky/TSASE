import numpy as np

def hypoboloid(chi_area, energy, alpha):

    scale_energy = energy / alpha    
    Un = np.sqrt(chi_area * chi_area + scale_energy * scale_energy - 1)

    return Un

def multipler(chi_area, energy, alpha):

    Un = (1+alpha)*chi_area*energy

    return Un
