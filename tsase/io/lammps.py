import numpy

import ase
from ase.data import *

def get_sym(tempmass):
    for i, d in enumerate(atomic_masses):   
        if tempmass == atomic_masses[i]:
            return i
        else:
            continue

def write_lammps(filename, atoms):
    f = open(filename, 'w')
    f.write(' lammps data file generated by tsase\n\n')
    f.write(' %d atoms\n\n' % len(atoms))
    symbols = atoms.get_chemical_symbols()
    elements = list(set(symbols))
    elements.sort()
    f.write(' %d atom types\n\n' % len(elements))
    f.write(' 0.0 %12.6f xlo xhi\n' % atoms.cell[0][0])
    f.write(' 0.0 %12.6f ylo yhi\n' % atoms.cell[1][1])
    f.write(' 0.0 %12.6f zlo zhi\n\n' % atoms.cell[2][2])
    f.write(' Masses\n\n')
    for i in range(len(elements)):
        f.write(' %d %12.6f\n' % (i+1, tsase.data.elements[elements[i]]['mass']))
    f.write('\n')
    f.write(' Atoms\n\n')
    for i in range(len(atoms)):
        typ = elements.index(symbols[i]) + 1
        r = atoms[i].position
        f.write(' %d %d %12.6f %12.6f %12.6f\n' % (i+1, typ, r[0], r[1], r[2]))
    f.write('\n')
    f.write(' Velocities\n\n')
    for i in range(len(atoms)):
        f.write(' %d 0.0 0.0 0.0\n' % (i+1))
    f.close()

def get_sym(tempmass):
    for i, d in enumerate(atomic_masses):
        if tempmass == atomic_masses[i]:
            return i
        else:
            continue

def read_lammps(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    index = 2
    while index < len(lines):
        if 'atoms' in lines[index]:
                natoms = int(lines[index].strip().split()[0])
                atoms = ase.Atoms(str(natoms)+'H')
                lenatomtypes = int(lines[index+1].strip().split()[0])
                index += 3               

        if 'xlo xhi' in lines[index]:
                atoms.cell[0][0] = float(lines[index].split()[1])-float(lines[index].split()[0])
                atoms.cell[1][1] = float(lines[index+1].split()[1])-float(lines[index+1].split()[0])                   
                atoms.cell[2][2] = float(lines[index+2].split()[1])-float(lines[index+2].split()[0])                    
                index += 4                    
        if 'Masses' in lines[index]:                    
                index +=2 
                idpoint={}
                for i in range(lenatomtypes):
                    typeid = int(lines[index].strip().split()[0])
                    tempmass = float(lines[index].strip().split()[1])
                    elnum = get_sym(tempmass)
                    idpoint[typeid] = elnum
                    index += 1
        index += 1
        if 'Atoms' in lines[index]:                    
                index +=2 
                for i in range(natoms):
                    data = lines[index].strip().split()                            
                    #yypa = atoms[int(data[0])-1]                            
                    idn = int(data[0])-1                            
                    idt = int(data[1])
                    atoms[idn].number = idpoint[idt]
                    atoms[idn].position[0] = float(data[2])                            
                    atoms[idn].position[1] = float(data[3])                            
                    atoms[idn].position[2] = float(data[4])                            
                    index += 1                            
    return atoms

def read_dump(filename):
    traj = []    
    lines = open(filename, 'r').readlines()
    index = 0
    while index < len(lines):
        if 'ITEM: TIMESTEP' in lines[index]:
            index += 2
        if 'ITEM: NUMBER OF ATOMS' in lines[index]:
            natoms = int(lines[index+1])
            atoms = ase.Atoms(str(natoms)+'H')
            traj.append(atoms)
            index += 2
        if 'ITEM: BOX BOUNDS' in lines[index]:
            atoms.cell[0][0] = float(lines[index+1].split()[1])
            atoms.cell[1][1] = float(lines[index+2].split()[1])
            atoms.cell[2][2] = float(lines[index+3].split()[1])
            index += 4
        if 'ITEM: ATOMS' in lines[index]:
            index += 1
            for i in range(natoms):
                data = lines[index].strip().split()
                a = atoms[int(data[0])-1]
                a.number = data[1]
                a.position[0] = float(data[2])
                a.position[1] = float(data[3])
                a.position[2] = float(data[4])
                index += 1
    return traj         
                
    
if __name__ == '__main__':
    import sys
    read_lammps(sys.argv[1])
#    traj = read_dump(sys.argv[1])
#    tsase.io.write_con(sys.argv[1] + '.con', traj[0], 'w')
#    for t in traj[1:]:
#        tsase.io.write_con(sys.argv[1] + '.con', t, 'a')
    
    
    
    
