#!/usr/bin/env python


import os
import sys
import numpy 
import glob
import shutil
import math

from optparse import OptionParser

from kdb import *

PBC_MAPPING_CHECK = False
REBOX_SUGGESTIONS = False
REMOVE_DUPLICATES = False
AUTOSCALE         = True


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
        c += a.r[i]
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
            sys.exit()
        a = working.pop()
        for i in undone[:]:
            v = pbc(temp.r[i] - temp.r[a], temp.box)
            d = numpy.linalg.norm(v)
            if d < (elements[temp.names[a]]["radius"] + elements[temp.names[i]]["radius"]) * (1.0 + NEIGHBOR_FUDGE):
                temp.r[i] = temp.r[a] + v
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
    elementPaths = [os.path.join(kdbdir, nameList) for nameList in namePermutations(sorted(reactant.getNameList()))]
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


if __name__ == "__main__":

    # Parse command line options.
    parser = OptionParser(usage = "%prog [options] reactant.con")
    parser.add_option("-d", "--kdbdir", dest = "kdbdir", 
                      help = "the path to the kinetic database",
                      default = "./kdb")
    parser.add_option("-a", "--af", dest = "af", action="store", type="float", 
                      help = "angstrom fudge parameter",
                      default = DISTANCE_CUTOFF)
    parser.add_option("-n", "--nf", dest = "nf", action="store", type="float", 
                      help = "neighbor fudge parameter",
                      default = NEIGHBOR_FUDGE)
    parser.add_option("--nodupes", dest = "nodupes", action="store_true",
                      help = "detect and remove duplicate suggestions (can be expensive)")
    options, args = parser.parse_args()
    DISTANCE_CUTOFF = options.af
    NEIGHBOR_FUDGE = options.nf
    REMOVE_DUPLICATES = options.nodupes

    # Make sure we get the reactant, saddle, product, and mode file names.
    if len(args) < 1:
        parser.print_help()
        sys.exit()
        
    # Load the reactant con file.
    reactant = loadcon(args[0])
    
    # Get the ibox to speed up pbcs.
    ibox = numpy.linalg.inv(reactant.box)
    
    # Create a directory for any suggestions we might make.
    if os.path.isdir("kdbmatches"):
        shutil.rmtree("kdbmatches")
    os.mkdir("kdbmatches")

    # A list of unique saddles, used for duplicate removal.
    uniques = []

    # Get a list of kdb entries that match the query configuration elementally.
    entries = getKDBentries(options.kdbdir, reactant)
    if len(entries) == 0:
        print "No entries for those elements."
        sys.exit()

    # For each nonfrozen atom in reactant, create a list of neighboring element
    # types and the count of each type.
    # TODO: this can be made N^2/2 trivially.
    # TODO: this can use SAP for ortho boxes.
    reactantNeighbors = {}
    for i in range(len(reactant)):
        if not reactant.free[i]:
            continue
        r1 = elements[reactant.names[i]]["radius"]
        reactantNeighbors[i] = {}
        for j in range(len(reactant)):
            if j == i:
                continue
            r2 = elements[reactant.names[j]]["radius"]
            d = numpy.linalg.norm(pbc(reactant.r[i] - reactant.r[j], reactant.box, ibox))
            if d > (r1 + r2) * (1 + NEIGHBOR_FUDGE):
                continue
            if reactant.names[j] not in reactantNeighbors[i]:
                reactantNeighbors[i][reactant.names[j]] = 0
            reactantNeighbors[i][reactant.names[j]] += 1
    
    # Create a list of element types and counts for the entire reactant. 
    reactantNameCount = reactant.nameCount()
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
        kdbmin = loadxyz(entry["minimum"])        

        # Make sure the reactant has at least as many atoms of each type as the
        # kdb configuration.
        passedNameCount = True
        kdbNameCount = kdbmin.nameCount()
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
                kdbmin.r[i] += 2.0 * (kdbmin.r[0] - kdbmin.r[i])
        
        # For each mobile atom in kdbmin, create a list of neighboring element
        # types and the count of each type.
        kdbNeighbors = {}
        for i in kdbmobile:
            r1 = elements[kdbmin.names[i]]["radius"]
            kdbNeighbors[i] = {}
            for j in range(len(kdbmin)):
                if j == i:
                    continue
                r2 = elements[kdbmin.names[j]]["radius"]
                d = numpy.linalg.norm(kdbmin.r[i] - kdbmin.r[j])
                if d > (r1 + r2) * (1 + NEIGHBOR_FUDGE):
                    continue
                if kdbmin.names[j] not in kdbNeighbors[i]:
                    kdbNeighbors[i][kdbmin.names[j]] = 0
                kdbNeighbors[i][kdbmin.names[j]] += 1

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
                kdbDistances[i] = numpy.linalg.norm(kdbmin.r[kdbAtom] - kdbmin.r[i])
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
                        pbcVector = reactant.atomAtomPbcVector(RA, reactantAtom)
                        if PBC_MAPPING_CHECK:
                            if not isDistance(pbcVector, kdbDistances[DA], reactant.box):
                                break
                        else:
                            if abs(kdbDistances[DA] - reactant.atomAtomPbcDistance(RA, reactantAtom)) > DISTANCE_CUTOFF:
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
            kdbrot.box = reactant.box.copy()
            
            # Rotation Matrix calculation start
            tb = kdbrot.copy()
            tb.r -= centroid(tb)
            ta = tb.copy()
            offset = centroid(reactantrot, mapping.values())
            i = 0
            for m in mapping:
                ta.r[i] = tb.r[m] + pbc((reactantrot.r[mapping[m]] - offset) - tb.r[m], reactantrot.box)
                i += 1
            ta.r -= centroid(ta)
            m = numpy.dot(tb.r.transpose(), ta.r)
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
            kdbrot.r -= translation1
            kdbrot.r = numpy.dot(kdbrot.r, Rmat)
            
            translation2 = centroid(reactantrot, mapping.values())
            
            kdbrot.r += translation2

            # Calculate a score for this mapping.
            score = max([numpy.linalg.norm(pbc(kdbrot.r[m] - reactantrot.r[mapping[m]], reactantrot.box)) for m in mapping])
            
            if score > DISTANCE_CUTOFF:
                continue
            
            # Load the saddle from the database.
            kdbSaddle = loadxyz(entry["saddle"])    
            
            # Mirror the saddle if the mirror flag is set for this entry.
            if entry["mirror"]:
                for i in range(len(kdbSaddle)):
                    kdbSaddle.r[i] += 2.0 * (kdbmin.r[0] - kdbSaddle.r[i])

            # Load the product from the database.
            kdbProduct = loadxyz(entry["product"])    
            
            # Mirror the product if the mirror flag is set for this entry.
            if entry["mirror"]:
                for i in range(len(kdbProduct)):
                    kdbProduct.r[i] += 2.0 * (kdbmin.r[0] - kdbProduct.r[i])

            # Map the mode.
            modeTemp = reactantrot.r * 0.0
            for m in mapping:
                modeTemp[mapping[m]] = mode[m]
            modeTemp /= numpy.linalg.norm(modeTemp)

            # Perform the saddle transformation.
            kdbSaddle.r -= translation1
            kdbSaddle.r = numpy.dot(kdbSaddle.r, Rmat)
            kdbSaddle.r += translation2
            
            # Perform the mode transformation.
            modeTemp = numpy.dot(modeTemp, Rmat)

            # Perform the product transformation.
            kdbProduct.r -= translation1
            kdbProduct.r = numpy.dot(kdbProduct.r, Rmat)
            kdbProduct.r += translation2
            
            # Create the suggestion.
            suggestion = reactant.copy()
            sugproduct = reactant.copy()
            for m in mapping:
                if suggestion.free[mapping[m]]:
                    suggestion.r[mapping[m]] = kdbSaddle.r[m]
                if sugproduct.free[mapping[m]]:
                    sugproduct.r[mapping[m]] = kdbProduct.r[m]
            
            # Check for duplicates.
            if REMOVE_DUPLICATES:
                isdupe = False
                for unique in uniques:
                    pan = per_atom_norm(unique.r - suggestion.r, suggestion.box, ibox)
                    if max(pan) <= DISTANCE_CUTOFF:
                        isdupe = True
                        break
                if isdupe:
                    continue
                uniques.append(suggestion.copy())     
            
            # Rebox.
            if REBOX_SUGGESTIONS:
                suggestion.r = pbc(suggestion.r, suggestion.box)
                sugproduct.r = pbc(sugproduct.r, sugproduct.box)
                
            # Write suggestion.
            savecon("kdbmatches/SADDLE_%d" % numMatches, suggestion)
            savecon("kdbmatches/PRODUCT_%d" % numMatches, sugproduct)
            save_mode("kdbmatches/MODE_%d" % numMatches, modeTemp)
            if os.path.isfile(entry["barrier"]): shutil.copyfile(entry["barrier"], "kdbmatches/BARRIER_%d" % numMatches)
            if os.path.isfile(entry["prefactor"]): shutil.copyfile(entry["prefactor"], "kdbmatches/PREFACTOR_%d" % numMatches)
            os.system("touch kdbmatches/.done_%d" % numMatches)                        
            
            #save debug xyz suggestion file.
            savexyz("kdbmatches/saddle_%d.xyz" % numMatches, suggestion)

            entryMatches += 1
            numMatches += 1

        print "%10d" % entryMatches
            

            
            





















