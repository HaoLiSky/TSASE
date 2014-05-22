
import numpy
import ase

def read_bopfox(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    line_index = 0
    newlines = []
    # strip out comment lines
    for line in lines:
        l = line.strip()
        if len(l) < 1:
          continue
        if l.startswith("/"):
          continue
        if l.startswith("#"):
          continue
        newlines.append(l)
    lines = newlines
    a = ase.Atoms()
    box = numpy.zeros((3,3))
    coord = "cartesian"
    scale = 1.0
    while line_index < len(lines):
        l = lines[line_index]
        if "magnetisation" in l:
            break
        if "=" in l:
            key = l.split("=")[0].strip()
            if key == "StrucName":
                pass
            elif key == "aLat":
                scale = float(l.split("=")[1].strip())
            elif key == "a1":
                box[0] = [float(i) for i in l.split("=")[1].strip().split()]
            elif key == "a2":
                box[1] = [float(i) for i in l.split("=")[1].strip().split()]
            elif key == "a3":
                box[2] = [float(i) for i in l.split("=")[1].strip().split()]
            elif key == "coord":
                coord = l.split("=")[1].strip()
        else:
            atomdata = l.strip().split()
            symbol = atomdata[0]
            x, y, z = float(atomdata[1]), float(atomdata[2]), float(atomdata[3])
            a.append(ase.Atom(symbol, (x, y, z)))
        line_index += 1
    a.cell = box * scale
    if coord == "direct":
        for i in range(len(a)):
            a.positions[i] = numpy.dot(a.positions[i], a.cell)
    else: 
        for i in range(len(a)):
            a.positions[i] *= scale
    return a

def write_bopfox(filename, a, w = 'w'):
    fout = open(filename, w)
    fout.write("StrucName = undefined\n")
    fout.write("aLat = 1.0\n")
    fout.write("a1 = % .12f   % .12f   % .12f\n" % (a.cell[0][0], a.cell[0][1], a.cell[0][2]))
    fout.write("a2 = % .12f   % .12f   % .12f\n" % (a.cell[1][0], a.cell[1][1], a.cell[1][2]))
    fout.write("a3 = % .12f   % .12f   % .12f\n" % (a.cell[2][0], a.cell[2][1], a.cell[2][2]))
    fout.write("coord = cartesian\n")
    symbols = a.get_chemical_symbols()
    r = a.get_positions()
    for i in range(len(a)):
        fout.write("%s % .12f % .12f % .12f /% 3d\n" % (symbols[i], r[i][0], r[i][1], r[i][2], i))
    fout.close()

