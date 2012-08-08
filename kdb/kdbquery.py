#!/usr/bin/env python


import os
import sys
import numpy 
import glob
import shutil
import math

from optparse import OptionParser

from kdb import *

from tsase.io import read_con, write_con
from ase.io.xyz import read_xyz, write_xyz

PBC_MAPPING_CHECK = False
REBOX_SUGGESTIONS = False
REMOVE_DUPLICATES = False

def isDistance(pbcvector, target, box):
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            for z in [-1, 0, 1]:
                temp = pbcvector.copy()
                temp += x * box[0]
                temp += y * box[1]
                temp += z * box[2]
                if abs(numpy.linalg.norm(temp) - target) < DISTANCE_CUTOFF:
                    return True
    return False
    

def centroid(a, which=None):
    if which == None:
        which = range(len(a))
    c = numpy.array([0.0, 0.0, 0.0])
    for i in which:
        c += a.positions[i]
    c /= len(which)
    return c


def clump(c, atoms):
    # Remove PBC's.
    temp = c.copy()
    undone = atoms[:]
    working = [undone.pop()]        
    while len(undone) > 0:
        if len(working) == 0:
            print "Dissociated reactant, or neighbor_fudge too small."
            return
        a = working.pop()
        for i in undone[:]:
            v = pbc(temp.positions[i] - temp.positions[a], temp.cell)
            d = numpy.linalg.norm(v)
            if d < (elements[temp.get_chemical_symbols()[a]]["radius"] + elements[temp.get_chemical_symbols()[i]]["radius"]) * (1.0 + NEIGHBOR_FUDGE):
                temp.positions[i] = temp.positions[a] + v
                working.append(i)
                undone.remove(i)
    return temp
    

def namePermutations(nameList):
    perms = [nameList[-1]]
    for i in range(len(nameList) - 2, -1, -1):
        for p in perms[:]:
            perms.append(nameList[i] + p)
        perms.append(nameList[i])
    return perms


def getKDBentries(kdbdir, reactant):
    elementPaths = [os.path.join(kdbdir, nameList) for nameList in namePermutations(sorted(getNameList(reactant)))]
    entries = []
    for elementPath in elementPaths:
        if not os.path.isdir(elementPath):
            continue
        for procDir in glob.glob(os.path.join(elementPath, "*")):
            entry = {"minimum": os.path.join(procDir, "min1.xyz"), 
                     "saddle": os.path.join(procDir, "saddle.xyz"), 
                     "mode": os.path.join(procDir, "mode"),
                     "mobile": os.path.join(procDir, "mobile"),
                     "product": os.path.join(procDir, "min2.xyz"),
                     "barrier": os.path.join(procDir, "barrier1"), 
                     "prefactor": os.path.join(procDir, "prefactor1"), 
                     "mirror": False}
            entries.append(entry)
            entry = {"minimum": os.path.join(procDir, "min1.xyz"), 
                     "saddle": os.path.join(procDir, "saddle.xyz"), 
                     "mode": os.path.join(procDir, "mode"),
                     "mobile": os.path.join(procDir, "mobile"),
                     "product": os.path.join(procDir, "min2.xyz"),
                     "barrier": os.path.join(procDir, "barrier1"), 
                     "prefactor": os.path.join(procDir, "prefactor1"), 
                     "mirror": True}
            entries.append(entry)
            entry = {"minimum": os.path.join(procDir, "min2.xyz"), 
                     "saddle": os.path.join(procDir, "saddle.xyz"), 
                     "mode": os.path.join(procDir, "mode"),
                     "mobile": os.path.join(procDir, "mobile"),
                     "product": os.path.join(procDir, "min1.xyz"),
                     "barrier": os.path.join(procDir, "barrier2"), 
                     "prefactor": os.path.join(procDir, "prefactor2"), 
                     "mirror": False}
            entries.append(entry)
            entry = {"minimum": os.path.join(procDir, "min2.xyz"), 
                     "saddle": os.path.join(procDir, "saddle.xyz"), 
                     "mode": os.path.join(procDir, "mode"),
                     "mobile": os.path.join(procDir, "mobile"),
                     "product": os.path.join(procDir, "min1.xyz"),
                     "barrier": os.path.join(procDir, "barrier2"), 
                     "prefactor": os.path.join(procDir, "prefactor2"), 
                     "mirror": True}
            entries.append(entry)
    return entries


def query(reactant, kdbdir, outputdir = "./kdbmatches", nf=0.2, dc=0.3, nodupes = False):
    global DISTANCE_CUTOFF, NEIGHBOR_FUDGE, REMOVE_DUPLICATES
    DISTANCE_CUTOFF = dc
    NEIGHBOR_FUDGE = nf
    REMOVE_DUPLICATES = nodupes

    # Get the ibox to speed up pbcs.
    ibox = numpy.linalg.inv(reactant.cell)
    
    # Create a directory for any suggestions we might make.
    if os.path.isdir(outputdir):
        shutil.rmtree(outputdir)
    os.mkdir(outputdir)

    # A list of unique saddles, used for duplicate removal.
    uniques = []

    # Get a list of kdb entries that match the query configuration elementally.
    entries = getKDBentries(kdbdir, reactant)
    if len(entries) == 0:
        print "No entries for those elements."
        return

    # For each nonfrozen atom in reactant, create a list of neighboring element
    # types and the count of each type.
    # TODO: this can be made N^2/2 trivially.
    # TODO: this can use SAP for ortho boxes.
    reactantNeighbors = {}
    for i in range(len(reactant)):
        if i in reactant.constraints[0].index:
            continue
        r1 = elements[reactant.get_chemical_symbols()[i]]["radius"]
        reactantNeighbors[i] = {}
        for j in range(len(reactant)):
            if j == i:
                continue
            r2 = elements[reactant.get_chemical_symbols()[j]]["radius"]
            d = numpy.linalg.norm(pbc(reactant.positions[i] - reactant.positions[j], reactant.cell, ibox))
            if d > (r1 + r2) * (1 + NEIGHBOR_FUDGE):
                continue
            if reactant.get_chemical_symbols()[j] not in reactantNeighbors[i]:
                reactantNeighbors[i][reactant.get_chemical_symbols()[j]] = 0
            reactantNeighbors[i][reactant.get_chemical_symbols()[j]] += 1
    
    # Create a list of element types and counts for the entire reactant. 
    reactantNameCount = nameCount(reactant)
    numMatches = 0
    
    ###########################################################################
    # (Main) Loop over each kdb entry.
    ###########################################################################
    for entry in entries:
    
        entryMatches = 0
    
        mirrored = "not mirrored"
        if entry["mirror"]:
            mirrored = "mirrored"
        print "checking %32s %16s" % (entry["minimum"], mirrored), 
        
        # Load the minimum.
        kdbmin = read_xyz(entry["minimum"])        

        # Make sure the reactant has at least as many atoms of each type as the
        # kdb configuration.
        passedNameCount = True
        kdbNameCount = nameCount(kdbmin)
        for name in kdbNameCount:
            if name not in reactantNameCount:
                passedNameCount = False
                break
            if kdbNameCount[name] > reactantNameCount[name]:
                passedNameCount = False
                break
        if not passedNameCount:
            print "%10d  name count fail" % entryMatches
            continue

        # Load the mobile atoms list.
        kdbmobile = [int(i.strip()) for i in open(entry["mobile"]).readlines()]        
        
        # Mirror the minimum if the mirror flag is set for this entry.
        if entry["mirror"]:
            for i in range(len(kdbmin)):
                kdbmin.positions[i] += 2.0 * (kdbmin.positions[0] - kdbmin.positions[i])
        
        # For each mobile atom in kdbmin, create a list of neighboring element
        # types and the count of each type.
        kdbNeighbors = {}
        for i in kdbmobile:
            r1 = elements[kdbmin.get_chemical_symbols()[i]]["radius"]
            kdbNeighbors[i] = {}
            for j in range(len(kdbmin)):
                if j == i:
                    continue
                r2 = elements[kdbmin.get_chemical_symbols()[j]]["radius"]
                d = numpy.linalg.norm(kdbmin.positions[i] - kdbmin.positions[j])
                if d > (r1 + r2) * (1 + NEIGHBOR_FUDGE):
                    continue
                if kdbmin.get_chemical_symbols()[j] not in kdbNeighbors[i]:
                    kdbNeighbors[i][kdbmin.get_chemical_symbols()[j]] = 0
                kdbNeighbors[i][kdbmin.get_chemical_symbols()[j]] += 1

        kdbUnmapped = range(len(kdbmin)) # Keep track of the kdb atoms that have been mapped.

        # Create the initial mappings.
        mappings = None
        db_a = kdbmobile[0] # This will be the selected mobile atom.
        for m in kdbmobile:
            mMappings = []
            for freeAtom in reactantNeighbors.keys():
                for elementType in reactantNeighbors[freeAtom]:
                    if elementType not in kdbNeighbors[m]:
                        break
                    if kdbNeighbors[m][elementType] != reactantNeighbors[freeAtom][elementType]:
                        break
                else:
                    mMappings.append({m:freeAtom})
            if mappings == None:
                mappings = mMappings
            if len(mMappings) < len(mappings):
                mappings = mMappings
                db_a = m
        
        kdbUnmapped.remove(db_a)
        
        while len(kdbUnmapped) > 0 and len(mappings) > 0:
            # Create a list of new mappings that will replace mappings at the
            # end of this iteration.
            newMappings = []
            # Select an unmapped atom from kdbmin.
            kdbAtom = kdbUnmapped.pop()
            # Get the distance between kdbAtom and every other atom in the kdb
            # configuration.
            kdbDistances = {}
            for i in range(len(kdbmin)):
                kdbDistances[i] = numpy.linalg.norm(kdbmin.positions[kdbAtom] - kdbmin.positions[i])
            # Loop over each mapping and try to place kdbAtom.
            for mapping in mappings:
                # Loop over each atom in the reactant.
                for reactantAtom in range(len(reactant)):
                    # Make sure it has not already been mapped.
                    if reactantAtom in mapping.values():
                        continue
                    # Loop over the atoms in mapping and see if the distance
                    # between reactantAtom and mapping.values() atoms is the same
                    # within DISTANCE_CUTOFF of the distance between kdbAtom
                    # and mapping.keys() atoms.
                    for DA in mapping.keys():
                        RA = mapping[DA]
                        pbcVector = atomAtomPbcVector(reactant, RA, reactantAtom)
                        if PBC_MAPPING_CHECK:
                            if not isDistance(pbcVector, kdbDistances[DA], reactant.cell):
                                break
                        else:
                            if abs(kdbDistances[DA] - atomAtomPbcDistance(reactant, RA, reactantAtom)) > DISTANCE_CUTOFF:
                                break
                    else:
                        newMapping = mapping.copy()
                        newMapping[kdbAtom] = reactantAtom
                        newMappings.append(newMapping)
            mappings = newMappings

        # Load the mode.
        mode = numpy.array([[float(item) for item in line.strip().split()] for line in open(entry["mode"], 'r').readlines()])

        # Loop over each mapping and try to find a rotation that aligns the
        # kdb configuration with the query configuration.
        for mapping in mappings:
        
            reactantrot = clump(reactant, mapping.values())
        
            # Make a copy of kdbmin for rotation and put it in the box.
            kdbrot = kdbmin.copy()
            kdbrot.cell = reactant.cell.copy()
            
            # Rotation Matrix calculation start
            tb = kdbrot.copy()
            tb.positions -= centroid(tb)
            ta = tb.copy()
            offset = centroid(reactantrot, mapping.values())
            i = 0
            for m in mapping:
                ta.positions[i] = tb.positions[m] + pbc((reactantrot.positions[mapping[m]] - offset) - tb.positions[m], reactantrot.cell)
                i += 1
            ta.positions -= centroid(ta)
            m = numpy.dot(tb.positions.transpose(), ta.positions)
            sxx = m[0][0]
            sxy = m[0][1]
            sxz = m[0][2]
            syx = m[1][0]
            syy = m[1][1]
            syz = m[1][2]
            szx = m[2][0]
            szy = m[2][1]
            szz = m[2][2]
            n = numpy.zeros((4,4))
            n[0][1] = syz - szy
            n[0][2] = szx - sxz
            n[0][3] = sxy - syx
            n[1][2] = sxy + syx
            n[1][3] = szx + sxz
            n[2][3] = syz + szy
            n += n.transpose()
            n[0][0] =  sxx + syy + szz
            n[1][1] =  sxx - syy - szz
            n[2][2] = -sxx + syy - szz
            n[3][3] = -sxx - syy + szz
            w, v = numpy.linalg.eig(n)
            maxw = 0
            maxv = 0
            for i in range(len(w)):
                if w[i] > maxw:
                    maxw = w[i]
                    maxv = v[:,i]
            Rmat = numpy.zeros((3,3))
            aa = maxv[0]**2
            bb = maxv[1]**2
            cc = maxv[2]**2
            dd = maxv[3]**2
            ab = maxv[0]*maxv[1]
            ac = maxv[0]*maxv[2]
            ad = maxv[0]*maxv[3]
            bc = maxv[1]*maxv[2]
            bd = maxv[1]*maxv[3]
            cd = maxv[2]*maxv[3]
            Rmat[0][0] = aa + bb - cc - dd
            Rmat[0][1] = 2*(bc-ad) 
            Rmat[0][2] = 2*(bd+ac) 
            Rmat[1][0] = 2*(bc+ad) 
            Rmat[1][1] = aa - bb + cc - dd
            Rmat[1][2] = 2*(cd-ab) 
            Rmat[2][0] = 2*(bd-ac) 
            Rmat[2][1] = 2*(cd+ab) 
            Rmat[2][2] = aa - bb - cc + dd
            Rmat = Rmat.transpose()
            # Rotation Matrix calculation end

            translation1 = centroid(kdbrot)
            kdbrot.positions -= translation1
            kdbrot.positions = numpy.dot(kdbrot.positions, Rmat)
            
            translation2 = centroid(reactantrot, mapping.values())
            
            kdbrot.positions += translation2

            # Calculate a score for this mapping.
            score = max([numpy.linalg.norm(pbc(kdbrot.positions[m] - reactantrot.positions[mapping[m]], reactantrot.cell)) for m in mapping])
            
            if score > DISTANCE_CUTOFF:
                continue
            
            # Load the saddle from the database.
            kdbSaddle = read_xyz(entry["saddle"])    
            
            # Mirror the saddle if the mirror flag is set for this entry.
            if entry["mirror"]:
                for i in range(len(kdbSaddle)):
                    kdbSaddle.positions[i] += 2.0 * (kdbmin.positions[0] - kdbSaddle.positions[i])

            # Load the product from the database.
            kdbProduct = read_xyz(entry["product"])    
            
            # Mirror the product if the mirror flag is set for this entry.
            if entry["mirror"]:
                for i in range(len(kdbProduct)):
                    kdbProduct.positions[i] += 2.0 * (kdbmin.positions[0] - kdbProduct.positions[i])

            # Map the mode.
            modeTemp = reactantrot.positions * 0.0
            for m in mapping:
                modeTemp[mapping[m]] = mode[m]
            modeTemp /= numpy.linalg.norm(modeTemp)

            # Perform the saddle transformation.
            kdbSaddle.positions -= translation1
            kdbSaddle.positions = numpy.dot(kdbSaddle.positions, Rmat)
            kdbSaddle.positions += translation2
            
            # Perform the mode transformation.
            modeTemp = numpy.dot(modeTemp, Rmat)

            # Perform the product transformation.
            kdbProduct.positions -= translation1
            kdbProduct.positions = numpy.dot(kdbProduct.positions, Rmat)
            kdbProduct.positions += translation2
            
            # Create the suggestion.
            suggestion = reactant.copy()
            sugproduct = reactant.copy()
            for m in mapping:
                if mapping[m] not in suggestion.constraints[0].index:
                    suggestion.positions[mapping[m]] = kdbSaddle.positions[m]
                if mapping[m] not in sugproduct.constraints[0].index:
                    sugproduct.positions[mapping[m]] = kdbProduct.positions[m]
            
            # Check for duplicates.
            if REMOVE_DUPLICATES:
                isdupe = False
                for unique in uniques:
                    pan = per_atom_norm(unique.positions - suggestion.positions, suggestion.cell, ibox)
                    if max(pan) <= DISTANCE_CUTOFF:
                        isdupe = True
                        break
                if isdupe:
                    continue
                uniques.append(suggestion.copy())     
            
            # Rebox.
            if REBOX_SUGGESTIONS:
                suggestion.positions = pbc(suggestion.positions, suggestion.cell)
                sugproduct.positions = pbc(sugproduct.positions, sugproduct.cell)
                
            # Write suggestion.
            write_con(outputdir + "/SADDLE_%d" % numMatches, suggestion)
            write_con(outputdir + "/PRODUCT_%d" % numMatches, sugproduct)
            save_mode(outputdir + "/MODE_%d" % numMatches, modeTemp)
            if os.path.isfile(entry["barrier"]): shutil.copyfile(entry["barrier"], outputdir + "/BARRIER_%d" % numMatches)
            if os.path.isfile(entry["prefactor"]): shutil.copyfile(entry["prefactor"], outputdir + "/PREFACTOR_%d" % numMatches)
            os.system("touch %s/.done_%d" % (outputdir, numMatches))                        
            
            #save debug xyz suggestion file.
            write_xyz(outputdir + "/saddle_%d.xyz" % numMatches, suggestion)

            entryMatches += 1
            numMatches += 1

        print "%10d" % entryMatches


if __name__ == "__main__":

    # Parse command line options.
    parser = OptionParser(usage = "%prog [options] reactant.con")
    parser.add_option("-d", "--kdbdir", dest = "kdbdir", 
                      help = "the path to the kinetic database",
                      default = "./kdb")
    parser.add_option("-c", "--dc", dest = "dc", action="store", type="float", 
                      help = "distance cutoff parameter",
                      default = DISTANCE_CUTOFF)
    parser.add_option("-n", "--nf", dest = "nf", action="store", type="float", 
                      help = "neighbor fudge parameter",
                      default = NEIGHBOR_FUDGE)
    parser.add_option("--nodupes", dest = "nodupes", action="store_true",
                      help = "detect and remove duplicate suggestions (can be expensive)")
    options, args = parser.parse_args()

    # Make sure we get the reactant file name.
    if len(args) < 1:
        parser.print_help()
        sys.exit()
        
    # Load the reactant con file.
    reactant = read_con(args[0])
    
    query(reactant, options.kdbdir, "./kdbmatches", options.dc, options.nf, options.nodupes)

            

            
            





















