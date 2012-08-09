#!/usr/bin/env python

import os
import sys
import numpy 
import glob

from optparse import OptionParser

from tsase.io import read_con

import ase
from ase.io.xyz import read_xyz, write_xyz

from kdb import *

def coordination_numbers(p, nf):
    nl = []
    for a in range(len(p)):
        nl.append([])
        for b in range(len(p)):
            if b != a:
                dist = numpy.linalg.norm(p.get_positions()[a] - p.get_positions()[b])        
                if dist < (elements[p.get_chemical_symbols()[a]]["radius"] + 
                           elements[p.get_chemical_symbols()[b]]["radius"]) * (1.0 + nf):
                    nl[a].append(b)
    return [len(l) for l in nl]

def getMappings(a, b, nf, dc, mappings = None):
    """ A recursive depth-first search for a complete set of mappings from atoms
        in configuration a to atoms in configuration b. Do not use the mappings
        argument, this is only used internally for recursion.
        
        Returns None if no mapping was found, or a dictionary mapping atom 
        indices a to atom indices b.
        
        Note: If a and b are mirror images, this function will still return a 
        mapping from a to b, even though it may not be possible to align them 
        through translation and rotation. """
    # If this is the top-level user call, create and loop through top-level
    # mappings.
    if mappings == None:
        # Find the least common coordination number in b.
        bCoordinations = coordination_numbers(b, nf)
        bCoordinationsCounts = {}
        for coordination in bCoordinations:
            if coordination in bCoordinationsCounts:
                bCoordinationsCounts[coordination] += 1
            else:
                bCoordinationsCounts[coordination] = 1
        bLeastCommonCoordination = bCoordinationsCounts.keys()[0]
        for coordination in bCoordinationsCounts.keys():
            if bCoordinationsCounts[coordination] < bCoordinationsCounts[bLeastCommonCoordination]:
                bLeastCommonCoordination = coordination
        # Find one atom in a with the least common coordination number in b. 
        # If it does not exist, return None.
        aCoordinations = coordination_numbers(a, nf)
        try:
            aAtom = aCoordinations.index(bLeastCommonCoordination)
        except ValueError:
            return None
        # Create a mapping from the atom chosen from a to each of the atoms with
        # the least common coordination number in b, and recurse.
        for i in range(len(bCoordinations)):
            if bCoordinations[i] == bLeastCommonCoordination:
                # Make sure the element types are the same.
                if a.get_chemical_symbols()[aAtom] != b.get_chemical_symbols()[i]:
                    continue
                mappings = getMappings(a, b, nf, dc, {aAtom:i})
                # If the result is not none, then we found a successful mapping.
                if mappings is not None:
                    return mappings
        # There were no mappings.        
        return None
    
    # This is a recursed invocation of this function.
    else:
        # Find an atom from a that has not yet been mapped.
        unmappedA = 0
        while unmappedA < len(a):
            if unmappedA not in mappings.keys():
                break
            unmappedA += 1
        # Calculate the distances from unmappedA to all mapped a atoms.
        distances = {}
        for i in mappings.keys():
            distances[i] = atomAtomDistance(a, unmappedA, i)
        
        # Loop over each unmapped b atom. Compare the distances between it and 
        # the mapped b atoms to the corresponding distances between unmappedA 
        # and the mapped atoms. If everything is similar, create a new mapping
        # and recurse.
        for bAtom in range(len(b)):
            if bAtom not in mappings.values():
                for aAtom in distances:
                    # Break if type check fails.
                    if b.get_chemical_symbols()[bAtom] != a.get_chemical_symbols()[unmappedA]:
                        break
                    # Break if distance check fails  
                    bDist = atomAtomDistance(b, bAtom, mappings[aAtom])
                    if abs(distances[aAtom] - bDist) > dc:
                        break
                else:
                    # All distances were good, so create a new mapping.
                    newMappings = mappings.copy()
                    newMappings[unmappedA] = bAtom
                    # If this is now a complete mapping from a to b, return it.
                    if len(newMappings) == len(a):
                        return newMappings
                    # Otherwise, recurse.
                    newMappings = getMappings(a, b, nf, dc, newMappings)
                    # Pass any successful mapping up the recursion chain. 
                    if newMappings is not None:
                        return newMappings     
        # There were no mappings.   
        return None 


def stripUnselectedAtoms(atoms, selected):
    """ Removes any atoms from atoms that are not in selected and returns a new
    structure and a mapping from atoms in the old structure to atoms in the new 
    structure. """
    src = atoms.copy()
    dest = atoms.copy()
    while len(dest) > 0:
        dest.pop()
    mapping = {}
    index = 0
    constraints = []
    for i in selected:
        mapping[i] = index
        index += 1
        if i in src.constraints[0].index: 
            constraints.append(index)
        dest.append(src[i])
    dest.set_constraint(ase.constraints.FixAtoms(constraints))
    return dest, mapping
    
    
def getProcessMobileAtoms(r, s, p, mac):
    """ Returns a list of atom indices that move more than mac 
    between reactant and saddle, saddle and product, or 
    reactant and product. If no atoms move more than mac, returns
    the atom that moves the most. """
    mobileAtoms = []
    reactant2saddle = per_atom_norm(s.positions - r.positions, s.get_cell())
    product2saddle = per_atom_norm(s.positions - p.positions, s.get_cell())
    reactant2product = per_atom_norm(p.positions - r.positions, s.get_cell())
    for i in range(len(s)):
        if max(reactant2saddle[i], product2saddle[i], reactant2product[i]) > mac:
            mobileAtoms.append(i)
    if len(mobileAtoms) == 0:
        mobileAtoms.append(list(reactant2product).index(max(reactant2product)))
    return mobileAtoms


def getProcessNeighbors(mobileAtoms, r, s, p, nf):
    """ Given a list mobile atoms, a reactant, saddle, and product, 
    returns a list of neighboring atoms according to the nf (NEIGHBOR_FUDGE)
    paramter."""
    neighborAtoms = []
    for atom in mobileAtoms:
        r1 = elements[s.get_chemical_symbols()[atom]]["radius"]
        for i in range(len(s)):
            if i in mobileAtoms or i in neighborAtoms:
                continue
            r2 = elements[s.get_chemical_symbols()[i]]["radius"]
            maxDist = (r1 + r2) * (1.0 + nf)
            if atomAtomPbcDistance(r, atom, i) < maxDist:
                neighborAtoms.append(i)
            elif atomAtomPbcDistance(s, atom, i) < maxDist:
                neighborAtoms.append(i)
            elif atomAtomPbcDistance(p, atom, i) < maxDist:
                neighborAtoms.append(i)
    return neighborAtoms


def insert(reactant, saddle, product, mode, kdbdir="./kdb", nf=0.2, dc=0.3, mac=0.7):

    mobileAtoms = getProcessMobileAtoms(reactant, saddle, product, mac)
        
    selectedAtoms = mobileAtoms + getProcessNeighbors(mobileAtoms, reactant, product, saddle, nf)
    
    # Quit if not enough selected atoms.
    if len(selectedAtoms) < 2:
        print "Too few atoms in process, or neighbor_fudge too small."
        return
            
    # Remove unselected atoms.
    reactant, mapping = stripUnselectedAtoms(reactant, selectedAtoms)
    saddle, mapping = stripUnselectedAtoms(saddle, selectedAtoms)
    product, mapping = stripUnselectedAtoms(product, selectedAtoms)

    # Update the mode.
    newMode = numpy.zeros((len(selectedAtoms), 3))
    for m in mapping:
        newMode[mapping[m]] = mode[m]
    mode = newMode

    # Remove PBC's.
    temp = reactant.copy()
    undone = range(len(temp))
    working = [undone.pop()]        
    while len(undone) > 0:
        if len(working) == 0:
            print "Dissociated reactant, or neighbor_fudge too small."
            return
        a = working.pop()
        for i in undone[:]:
            v = pbc(temp.positions[i] - temp.positions[a], temp.get_cell())
            d = numpy.linalg.norm(v)
            if d < (elements[temp.get_chemical_symbols()[a]]["radius"] + 
                    elements[temp.get_chemical_symbols()[i]]["radius"]) * (1.0 + nf):
                temp[i].position = temp[a].position + v
                working.append(i)
                undone.remove(i)
    v1s = pbc(saddle.positions - reactant.positions, reactant.get_cell())
    v12 = pbc(product.positions - reactant.positions, reactant.get_cell())
    reactant = temp
    saddle.positions = reactant.positions + v1s
    product.positions = reactant.positions + v12
    
    # Find saddle center of coordinates.
    coc = numpy.zeros((1,3))
    for i in range(len(saddle)):
        coc += saddle[i].position
    coc = coc / len(saddle)
    
    # Shift all structures so that the saddle center of coordinates is at 
    # [0, 0, 0].
    reactant.positions = reactant.positions - coc    
    saddle.positions = saddle.positions - coc    
    product.positions = product.positions - coc    

    # Give all structures a huge box.
    # TODO: all references to boxes should be removed after PBCs are removed.
    reactant.cell = numpy.identity(3) * 1024
    saddle.cell = numpy.identity(3) * 1024
    product.cell = numpy.identity(3) * 1024

    # Get the element path for this process
    elementPath = "".join(getNameList(reactant))
    elementPath = os.path.join(kdbdir, elementPath)
    if not os.path.exists(elementPath):
        os.makedirs(elementPath)

    # Get a list of process subdirectories in the element path.
    procdirs = glob.glob(os.path.join(elementPath, "*"))

    # Loop over the existing process and check for matches to the saddle.
    for procdir in procdirs:
        dbSaddle = read_xyz(os.path.join(procdir, "saddle.xyz"))
        if len(saddle) != len(dbSaddle):
            continue
        if getMappings(saddle, dbSaddle, nf, dc) is not None:
            print "duplicate of", procdir
            return      

    # Create the path for this process.
    i = 0
    while os.path.exists(os.path.join(elementPath, str(i))):
        i += 1
    processPath = os.path.join(elementPath, str(i))
    os.makedirs(processPath)

    # Save the configurations for this process.  
    write_xyz(os.path.join(processPath, "min1.xyz"), reactant)
    write_xyz(os.path.join(processPath, "saddle.xyz"), saddle)
    write_xyz(os.path.join(processPath, "min2.xyz"), product)
    save_mode(os.path.join(processPath, "mode"), mode)

    def numberfile(filename, number):
        f = open(filename, 'w')
        f.write(str(number))
        f.close()
        
    # Save the list of mobile atoms.
    f = open(os.path.join(processPath, "mobile"), 'w')
    for atom in mobileAtoms:
        f.write("%d\n" % mapping[atom])
    f.close()
    
    # Save a movie of the local process.
    steps = 8
    movie = [reactant.copy()]
    for i in range(1, steps):
        temp = reactant.copy()
        temp.positions = reactant.positions + (saddle.positions - reactant.positions) * (i / float(steps))
        movie.append(temp)
    movie.append(saddle.copy())
    for i in range(1, steps):
        temp = reactant.copy()
        temp.positions = saddle.positions + (product.positions - saddle.positions) * (i / float(steps))
        movie.append(temp)
    movie.append(product.copy())
    write_xyz(os.path.join(processPath, "movie.xyz"), movie)
        
    # Indicate that the process was inserted successfully.
    print "good"


if __name__ == "__main__":

    # Parse command line options.
    parser = OptionParser(usage = "%prog [options] reactant.con saddle.con product.con mode")
    parser.add_option("-d", "--kdbdir", dest = "kdbdir", 
                      help = "the path to the kinetic database",
                      default = "./kdb")
    parser.add_option("-n", "--nf", dest = "nf", action="store", type="float", 
                      help = "neighbor fudge parameter",
                      default = NEIGHBOR_FUDGE)
    parser.add_option("-c", "--dc", dest = "dc", action="store", type="float", 
                      help = "distance cutoff parameter",
                      default = DISTANCE_CUTOFF)
    parser.add_option("-m", "--mac", dest = "mac", action="store", type="float", 
                      help = "mobile atom cutoff parameter",
                      default = MOBILE_ATOM_CUTOFF)
    options, args = parser.parse_args()

    # Make sure we get the reactant, saddle, product, and mode files.
    if len(args) < 4:
        parser.print_help()
        sys.exit()

    # Load the reactant, saddle, product, and mode files.
    reactant = read_con(args[0])
    saddle = read_con(args[1])
    product = read_con(args[2])
    mode = load_mode(args[3])
    
    insert(reactant, saddle, product, mode, options.kdbdir, 
           options.nf, options.dc, options.mac)

    
        
    
    
                    
                
                
        
        
    
    
    
        
        
        
        
        
            

        











