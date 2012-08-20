
import numpy

MOBILE_ATOM_CUTOFF = 0.7
NEIGHBOR_FUDGE = 0.2
DISTANCE_CUTOFF = 0.3


def atomAtomPbcVector(atoms, a, b):
    if not hasattr(atoms, 'ibox'):
        atoms.ibox = numpy.linalg.inv(atoms.get_cell())
    if not hasattr(atoms, 'pbcVectors'):
        atoms.pbcVectors = {}
    if (a, b) not in atoms.pbcVectors or (b, a) not in atoms.pbcVectors:
        atoms.pbcVectors[(a, b)] = pbc(atoms.positions[b] - atoms.positions[a], atoms.get_cell(), atoms.ibox)
        atoms.pbcVectors[(b, a)] = -atoms.pbcVectors[(a, b)]
    return atoms.pbcVectors[(a, b)]
            

def atomAtomPbcDistance(atoms, a, b):
    if not hasattr(atoms, 'pbcDistances'):
        atoms.pbcDistances = {}
    if (a, b) not in atoms.pbcDistances or (b, a) not in atoms.pbcDistances:
        atoms.pbcDistances[(a, b)] = numpy.linalg.norm(atomAtomPbcVector(atoms, a, b))
        atoms.pbcDistances[(b, a)] = atoms.pbcDistances[(a, b)]
    return atoms.pbcDistances[(a, b)]
        

def atomAtomDistance(atoms, a, b):
    if not hasattr(atoms, 'distances'):
        atoms.distances = {}
    if (a, b) not in atoms.distances or (b, a) not in atoms.distances:
        atoms.distances[(a, b)] = numpy.linalg.norm(atoms.positions[a] - atoms.positions[b])
        atoms.distances[(b, a)] = atoms.distances[(a, b)]
    return atoms.distances[(a, b)]


def getNameList(atoms):
    """
    Returns a sorted list of element names.
    """
    nl = []
    for name in atoms.get_chemical_symbols():
        if name not in nl:
            nl.append(name)
    return sorted(nl)


def nameCount(atoms):
    counts = {}
    for name in atoms.get_chemical_symbols():
        if not name in counts:
            counts[name] = 0
        counts[name] += 1
    return counts
        

def pbc(r, box, ibox = None):
    """
    Applies periodic boundary conditions.
    Parameters:
        r:      the vector the boundary conditions are applied to
        box:    the box that defines the boundary conditions
        ibox:   the inverse of the box. This will be calcluated if not provided.
    """
    if ibox == None:    
        ibox = numpy.linalg.inv(box)
    vdir = numpy.dot(r, ibox)
    vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
    return numpy.dot(vdir, box)


def per_atom_norm(v, box, ibox = None):
    '''
    Returns a length N numpy array containing per atom distance
        v:      an Nx3 numpy array
        box:    box matrix that defines the boundary conditions
        ibox:   the inverse of the box. will be calculated if not provided
    '''
    diff = pbc(v, box, ibox)
    return numpy.array([numpy.linalg.norm(d) for d in diff])


def load_mode(modefilein):
    ''' 
    Reads a mode.dat file into an N by 3 numpy array
        modefilein: filename
    '''
    f = open(modefilein, 'r')
    lines = f.readlines()
    f.close()
    mode = []
    for line in lines:
        l = line.strip().split()
        for j in range(3):
            mode.append(float(l[j]))
    mode = numpy.array(mode)
    mode.resize(len(mode)/3, 3)
    return mode


def save_mode(modefileout, displace_vector):
    '''
    Saves an Nx3 numpy array into a mode.dat file. 
        modefileout:     filename
        displace_vector: the mode (Nx3 numpy array)
    '''
    f = open(modefileout, 'w')
    for i in range(len(displace_vector)):
        f.write("%.3f %.3f %.3f\n" % (displace_vector[i][0], 
            displace_vector[i][1], displace_vector[i][2]))