##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

'''
Con(figuration) i/o library
'''
import ConfigParser
from cStringIO import StringIO
import logging
logger = logging.getLogger('io')
import numpy
import os

import cPickle as pickle

import atoms
import config

def save_prng_state():
    state = numpy.random.get_state()
    fh = open('prng.pkl', 'wb')
    pickle.dump(state, fh, pickle.HIGHEST_PROTOCOL)

def get_prng_state():
    fh = open('prng.pkl')
    state = pickle.load(fh)
    numpy.random.set_state(state)

def length_angle_to_box(boxlengths, angles):
    box = numpy.zeros( (3,3) )
    angles *= numpy.pi/180.0
    box[0][0] = 1.0
    box[1][0] = numpy.cos(angles[0])
    box[1][1] = numpy.sin(angles[0])
    box[2][0] = numpy.cos(angles[1])
    box[2][1] = (numpy.cos(angles[2])-box[1][0]*box[2][0])/box[1][1]
    box[2][2] = numpy.sqrt(1.0-box[2][0]**2-box[2][1]**2)
    box[0,:]*=boxlengths[0]
    box[1,:]*=boxlengths[1]
    box[2,:]*=boxlengths[2]
    return box

def box_to_length_angle(box):
    lengths = numpy.zeros(3)
    lengths[0] = numpy.linalg.norm(box[0,:])
    lengths[1] = numpy.linalg.norm(box[1,:])
    lengths[2] = numpy.linalg.norm(box[2,:])
    angles = numpy.zeros(3)
    angles[0] = numpy.arccos(numpy.dot(box[0,:]/lengths[0],box[1,:]/lengths[1]))
    angles[1] = numpy.arccos(numpy.dot(box[0,:]/lengths[0],box[2,:]/lengths[2]))
    angles[2] = numpy.arccos(numpy.dot(box[1,:]/lengths[1],box[2,:]/lengths[2]))
    angles *= 180.0/numpy.pi
    return lengths, angles


def loadcons(filename):
    filein = open(filename, 'r')
    p = []
    while True:
        try:
            p.append(loadcon(filein, reset=False))
        except:
            return p


def loadposcars(filename):
    filein = open(filename, 'r')
    p = []
    while True:
        try:
            p.append(loadposcar(filein))
        except:
            return p


def loadcon(filein, reset = True):
    '''
    Load a con file
        filein: may be either a filename or a file-like object
    '''
    if hasattr(filein, 'readline'):
        con = filein
    else:
        con = open(filein, 'r')
    con.readline() # line 1: comment
    con.readline() # line 2: comment
    # determine how many dimensions
    tmp = numpy.array(con.readline().split()) # line 3: Box lengths
    for i in range(len(tmp)):
        dim=i+1
        try: float(tmp[i])
        except:
            dim=i
            break
    # handle the box
    boxlengths=numpy.zeros(dim)
    for i in range(dim):
        boxlengths[i]=float(tmp[i])
    boxangles=numpy.array([ float(f) for f in con.readline().split()[0:dim] ]) # line 4: Box angles
    boxtemp=numpy.zeros((dim,dim),'d')
    boxtemp = length_angle_to_box(boxlengths,boxangles)
    con.readline() # line 5: comment
    con.readline() # line 6: comment
    num_types = int(con.readline().split()[0]) # line 7: number of atom types
    num_each_type = con.readline().split() # line 8: number of each type of atom
    mass_of_type = con.readline().split() # line 9: mass of each type of atom
    num_atoms = 0
    for i in range(num_types):
        num_each_type[i] = int(num_each_type[i])
        mass_of_type[i] = float(mass_of_type[i])
        num_atoms += num_each_type[i]
    a = atoms.Atoms(num_atoms)
    a.box = boxtemp
    index = 0
    for i in range(num_types):
        name = con.readline().strip()
        if abs(1.0-mass_of_type[i]) < 1e-6 and name != "H":
            logger.warning("WARNING: Mass of %s set to 1.0", name)

        con.readline() # skip meaningless line
        for j in range(num_each_type[i]):
            vals = con.readline().split()
            for k in range(dim):
                a.r[index][k] = float(vals[k])
            a.mass[index] = mass_of_type[i]
            a.names[index] = name
            if not int(vals[dim])==0:
                a.free[index]=0
            index += 1
    if reset:
        con.seek(0)
    return a

def savecon(fileout, p, w = 'w'):
    '''
    Save a con file
        fileout: can be either a file name or a file-like object
        p:       information (in the form of an atoms object) to save
        w:       write/append flag
    '''
    if hasattr(fileout, 'write'):
        con = fileout
    else:
        con = open(fileout, w)
    print >> con, "Generated by eOn"
    print >> con
    dim = len(p.r[0])
    lengths, angles = box_to_length_angle(p.box)
    print >> con, " ".join(['%.6f' % s for s in lengths])
    print >> con, " ".join(['%.6f' % s for s in angles])
    print >> con
    print >> con
    atom_count = {}
    name_order = []
    for i in range(len(p)):
        name = p.names[i]
        if name not in name_order:
            name_order.append(name)
        if name in atom_count:
            atom_count[name] += 1
        else:
            atom_count[name] = 1
    print >> con, len(name_order)
    print >> con, " ".join([str(atom_count[i]) for i in name_order])
    printmasses = []
    index = 0
    for i in range(len(name_order)):
        printmasses.append(p.mass[index])
        index += atom_count[name_order[i]]
    print >> con, " ".join(["%.6f"% i for i in printmasses])
    index = 0
    for i in range(len(name_order)):
        print >> con, name_order[i]
        print >> con, "Coordinates of component", i+1
        for j in range(atom_count[name_order[i]]):
            con.write("%.6f %.6f %.6f %d %d\n" %( p.r[index][0], p.r[index][1], p.r[index][2], int(not p.free[index]), index+1))
            index += 1


def load_mode(modefilein):
    ''' 
    Reads a mode.dat file into an N by 3 numpy array
        modefilein: may be either a file-like object of a filename
    '''
    if hasattr(modefilein, 'readline'):
        f = modefilein
    else:
        f = open(modefilein, 'r')
    if len(f.readline().split()) == 3:
        f.seek(0);
    lines = f.readlines()
    mode = []
    for line in lines:
        l = line.strip().split()
        for j in range(3):
            mode.append(float(l[j]))
    return numpy.array(mode).reshape(len(mode)/3, 3)

def save_mode(modefileout, displace_vector):
    '''
    Saves an Nx3 numpy array into a mode.dat file.
        modefileout:     may be either a filename or file-like object
        displace_vector: the mode (Nx3 numpy array)
    '''
    if hasattr(modefileout, 'write'):
        f = modefileout
    else:
        f = open(modefileout, 'w')
    for i in range(len(displace_vector)):
        f.write("%.3f %.3f %.3f\n" % (displace_vector[i][0], 
            displace_vector[i][1], displace_vector[i][2]))


def save_results_dat(fileout, results):
    '''
    Saves a results.dat file from a dictionary
    '''
    if hasattr(fileout, 'write'):
        f = fileout
    else:
        f = open(fileout, 'w')

    for key in results:
        print >> f, results[key], key

def modify_config(config_path, changes):
    parser = ConfigParser.SafeConfigParser()
    parser.read(config.config_path)
    for change in changes:
        parser.set(*change)
    config_str_io = StringIO()
    parser.write(config_str_io)
    config_str_io.seek(0)
    return config_str_io

def parse_results(filein):
    '''
    Reads a results.dat file and gives a dictionary of the values contained therein
    '''
    if hasattr(filein, 'readline'):
        f = filein
        f.seek(0)
    else:
        f = open(filein)
    results = {}
    for line in f:
        line = line.split()
        if len(line) < 2:
            continue
        if '.' in line[0]:
            try:
                results[line[1]] = float(line[0])
            except ValueError:
                logger.warning("Couldn't parse float in results.dat: %s", line)
        else:
            try:
                results[line[1]] = int(line[0])
            except ValueError:
                try:
                    results[line[1]] = line[0]
                except ValueError:
                    logger.warning("Couldn't parse string in results.dat: %s", line)

    return results

def loadposcar(filein):
    '''
    Load the POSCAR file named filename and returns an atoms object
    '''
    if hasattr(filein, 'readline'):
        f = filein
    else:
        f = open(filein, 'r')
    # Line 1: Atom types
    AtomTypes = f.readline().split() 
    # Line 2: scaling of coordinates
    scale = float(f.readline()) 
    # Lines 3-5: the box
    box = numpy.zeros((3, 3))
    for i in range(3):
        line = f.readline().split()
        box[i] = numpy.array([float(line[0]), float(line[1]), float(line[2])]) * scale
    # Line 6: number of atoms of each type. 
    line = f.readline().split()
    NumAtomsPerType = []
    for l in line:
        NumAtomsPerType.append(int(l))
    # Now have enough info to make the atoms object.
    num_atoms = sum(NumAtomsPerType)
    p = atoms.Atoms(num_atoms)
    # Fill in the box.
    p.box = box
    # Line 7: selective or cartesian
    sel = f.readline()[0]
    selective_flag = (sel == 's' or sel == 'S')
    if not selective_flag: 
        car = sel
    else:
        car = f.readline()[0]
    direct_flag = not (car == 'c' or car == 'C' or car == 'k' or car == 'K')
    atom_index = 0
    for i in range(len(NumAtomsPerType)):
        for j in range(NumAtomsPerType[i]):
            p.names[atom_index] = AtomTypes[i]
            line = f.readline().split()
            if(selective_flag):
                assert len(line) >= 6
            else:
                assert len(line) >= 3
            pos = line[0:3]
            if selective_flag:
                sel = line[3:7]
                if sel[0] == 'T' or sel[0] == 't': 
                    p.free[atom_index] = 1
                elif sel[0] == 'F' or sel[0] == 'f':
                    p.free[atom_index] = 0
            p.r[atom_index] = numpy.array([float(q) for q in pos])
            if direct_flag:
                p.r[atom_index] = numpy.dot(p.r[atom_index], p.box)
            else:
                p.r[atom_index] *= scale
            atom_index += 1
    return p


def saveposcar(fileout, p, w='w', direct = False):
    '''
    Save a POSCAR
        fileout: name to save it under
        point:    atoms object to save
        w:        write/append flag
    ''' 
    if hasattr(fileout, 'write'):
        poscar = fileout
    else:
        poscar = open(fileout, w)
    atom_types = []
    num_each_type = {}
    for name in p.names:
        if not name in atom_types:
            atom_types.append(name)
            num_each_type[name] = 1
        else:
            num_each_type[name] += 1
    print >> poscar, " ".join(atom_types) #Line 1: Atom types
    print >> poscar, 1.0 #Line 2: scaling
    for i in range(3):
        print >> poscar, " ".join(['%20.14f' % s for s in p.box[i]]) #lines 3-5: box
    print >> poscar, " ".join(['%s' % num_each_type[key] for key in atom_types])
    print >> poscar, 'Selective Dynamics' #line 6: selective dynamics
    if direct:
        print >> poscar, 'Direct' #line 7 cartesian coordinates
        ibox = numpy.linalg.inv(numpy.array(p.box))
        p.r = numpy.dot(p.r, ibox)
    else:
        print >> poscar, 'Cartesian' #line 7 cartesian coordinates
    for i in range(len(p)):
            posline = " ".join(['%20.14f' % s for s in p.r[i]]) + " "
            for j in range(3):
                if(p.free[i]): 
                    posline+='   T'
                else: 
                    posline+='   F'
            print >> poscar, posline


from ConfigParser import SafeConfigParser as SCP
class ini(SCP):

    def __init__(self, filenames):
        self.loaded = False
        self.filenames = filenames
        SCP.__init__(self)

    def read(self):
        self.loaded = True
        SCP.read(self, self.filenames)

    def get(self, section, option, default="ini_no_default"):
        if not self.loaded:
            self.read()
        try:
            value = SCP.get(self, section, option)
        except:
            if default == "ini_no_default":
                raise NameError("Section or option missing, no default specified")
            return default
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        if value.lower() == 'true':
            return True
        if value.lower() == 'false':
            return False
        return value

    def getint(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")
    def getfloat(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")
    def getboolean(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")

    def set(self, section, option, value):
        if section not in self.sections():
            self.add_section(section)
        SCP.set(self, section, option, str(value))
        if type(self.filenames) == str:  
            name = self.filenames
        else:
            name = self.filenames[-1]
        configfile = open(name, 'wb')
        self.write(configfile)
        configfile.close()


class Dynamics:
    """ The Dynamics class handles I/O for the dynamics.txt file of an aKMC simulation. """

    def __init__(self, filename):
        self.filename = filename
        if not os.path.exists(filename):
            f = open(self.filename, 'w')
            header = "%12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n" % ('step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate')
            f.write(header)
            f.write("-" * len(header))
            f.write("\n")
            f.close()
            self.next_step = 0

        # read last lines of the file to determine iteration nr
        else:
            f = open(self.filename,'r')
            f.seek(0,2)	#seek to EOF
            fsize = f.tell()
            # seek 1024 bytes back (or to beginning of file if fsize < 1024 )
            # last line must be contained in this block
            f.seek( max( fsize - 1024 , 0 ) , 0)
            lines = f.readlines()
            self.next_step = int ( lines[-1].split()[0] ) + 1 # determine iteration nr of next step

    def append(self, reactant_id, process_id, product_id, step_time, total_time, barrier, rate):
        f = open(self.filename, 'a')
        f.write("%12d  %12d  %12d  %12d  %12e  %12e  %12f  %12e\n" % (self.next_step, reactant_id, process_id, product_id, step_time, total_time, barrier, rate))
        f.close()
        self.next_step += 1

    def append_sb(self, reactant_id, process_id, product_id, step_time, total_time, basin_id):
        f = open(self.filename, 'a')
        f.write("%12d  %12d  %12d  %12d  %12e  %12e  %12d\n" % (self.next_step, reactant_id, process_id, product_id, step_time, total_time, basin_id))
        f.close()
        self.next_step += 1

    def get(self):
        f = open(self.filename, 'r')
        lines = f.readlines()[2:]
        f.close()
        data = []
        for line in lines:
            split = line.split()
            data.append({"reactant":    int(split[1]),
                         "process":     int(split[2]),
                         "product":     int(split[3]),
                         "steptime":    float(split[4]),
                         "totaltime":   float(split[5]),
                         "barrier":     float(split[6]),
                         "prefactor":   float(split[7])})
        return data

def load_potfiles(pot_dir):
    ret = {}
    if os.path.isdir(pot_dir):
        for name in os.listdir(pot_dir):
            if os.path.isdir(name):
                continue
            a = open(os.path.join(pot_dir, name), 'r')
            b = StringIO("".join(a.readlines()))
            ret[name] = b
    return ret

class TableException(Exception):
    pass

class Table:
    """
    A class that provides a nice io abstraction for table like data.  The data
    is saved in a pretty printed format. Also provides nice data retrival
    methods.

    >>> t = Table("sample.tbl", ['id', 'name', 'age' ])
    >>> t.eagerwrite = False
    >>> t.add_row({'id':0,'name':"Sam","age":24})
    >>> t.add_row({'id':1,'name':"David","age":50})
    >>> t.add_row({'id':2,'name':"Anna","age":21})
    >>> t #doctest: +NORMALIZE_WHITESPACE
        id name  age
        -- ----- ---
        0  Sam   24
        1  David 50
        2  Anna  21

    Rows can be accessed directly:
    >>> t.rows[1] #doctest: +SKIP
        {'age': 50, 'id': 1, 'name': 'David'}
    >>> t.max_value('age') #doctest: +NORMALIZE_WHITESPACE
        50
    >>> t.min_row('age') #doctest: +NORMALIZE_WHITESPACE +SKIP
        {'age': 21, 'id': 2, 'name': 'Anna'}
    >>> sorted(t.min_row('id').items()) #doctest: +NORMALIZE_WHITESPACE
        [('age', 24), ('id', 0), ('name', 'Sam')]
    >>> len(t) #doctest: +NORMALIZE_WHITESPACE
        3
    >>> sum(t.getcolumn('age')) #doctest: +NORMALIZE_WHITESPACE
        95
    >>> t.write() #doctest: +SKIP

    The table can be loaded from disk without specifying columns. This is
    slightly unsafe because the columns can't be checked, but it could cut down
    on the verbosity in some places.
    >>> t2 = Table("sample.tbl") #doctest: +SKIP
"""

    #XXX: This is the number of digits that a floating point number gets
    #     serialized with. Should it be some sort of config option?
    #     Or is there just a good default?

    def __init__(self, filename, columns=None, overwrite=False):
        self.filename = filename
        self.columns = columns
        self.rows = []
        self.columntypes = {}
        self.columnwidths = {}
        self.initialized = False
        self.overwrite = overwrite

        self.floatprecision = 6
        self.eagerwrite = True

    def init(self):
        """Checks to see if self.filename exists. If it does self.rows
           will be initialized from disk."""
        self.initialized = True
        if os.path.isfile(self.filename) and not self.overwrite:
            self.read(self.filename)
        else:
            if self.columns is None:
                raise TableException("Columns are not optional for new tables")

            for c in self.columns:
                self.columnwidths[c] = len(c)

    def read(self, filename):
        self.eagerwrite = False
        f = open(self.filename, "r")
        filecolumns = f.readline().split()
        if self.columns != None:
            if filecolumns != self.columns:
                raise TableException("Column name mismatch: %s" % filename)
        else:
            self.columns = filecolumns

        for c in self.columns:
            self.columnwidths[c] = len(c)

        # skip comment line
        f.readline()

        for line in f:
            fields = line.split()
            row = {}
            coli = 0
            for field in fields:
                try:
                    field = int(field)
                except ValueError:
                    try:
                        field = float(field)
                    except ValueError:
                        field = field.strip()
                        pass
                row[self.columns[coli]] = field
                coli += 1
            self.add_row(row)
        f.close()
        self.eagerwrite = True

    def __repr__(self):
        if not self.initialized:
            self.init()
        f = StringIO()
        self.writefilehandle(f)
        return f.getvalue()

    def __len__(self):
        if not self.initialized:
            self.init()
        return len(self.rows)

    def __iter__(self):
        for row in self.rows:
            yield row

    def write(self):
        if not self.initialized:
            self.init()
        f = open(self.filename, "w")
        self.writefilehandle(f)
        f.close()

    def writefilehandle(self, filehandle):
        f = filehandle
        line = ' '.join([ "%-*s"%(self.columnwidths[c], c) for c in self.columns ])
        f.write(line+"\n")

        line = ''
        for c in self.columns:
            line += '-'*self.columnwidths[c]+' '
        f.write(line+'\n')

        for row in self.rows:
            line = ""
            for c in self.columns:
                if self.columntypes[c] == float:
                    line += "%#-*.*G " % (self.columnwidths[c],self.floatprecision,
                                         row[c])
                else:
                    line += "%-*s " % (self.columnwidths[c],row[c])
            f.write(line+"\n")

    def add_row(self, row):
        if not self.initialized:
            self.init()
        mismatched_columns = set(self.columns).symmetric_difference(set(row.keys()))
        if len(mismatched_columns) != 0:
            raise TableException("Mismatched columns %s" % str(mismatched_columns))

        if len(self.rows) == 0:
            for c in row:
                self.columntypes[c] = type(row[c])
        else:
            for c in row:
                if type(row[c]) != self.columntypes[c]:
                    raise TableException("Type mismatch for column %s" % c)

        for c in row:
            if self.columntypes[c] == float:
                self.columnwidths[c] = max(self.columnwidths[c], self.floatprecision+5)
            else:
                self.columnwidths[c] = max(self.columnwidths[c], len(str(row[c])))

        self.rows.append(row)
        if self.eagerwrite:
            self.write()

    def delete_row(self, column, value):
        if not self.initialized:
            self.init()
        rows_to_delete = []
        for row in self.rows:
            if row[column] == value:
                rows_to_delete.append(row)
        map(self.rows.remove, rows_to_delete)
        if self.eagerwrite:
            self.write()
        return len(rows_to_delete)

    def delete_row_func(self, column, func):
        if not self.initialized:
            self.init()

        rows_to_delete = []
        for row in self.rows:
            if func(row[column]):
                rows_to_delete.append(row)
        map(self.rows.remove, rows_to_delete)
        if self.eagerwrite:
            self.write()
        return len(rows_to_delete)

    def find_value(self, column, func):
        if not self.initialized:
            self.init()
        value = None
        for row in self.rows:
            if value is None:
                value = row[column]
                continue
            value = func(value, row[column])
        return value

    def find_row(self, column, func):
        if not self.initialized:
            self.init()
        value = None
        for row in self.rows:
            if value is None:
                value = row
                continue

            if func(row[column],value[column])==row[column]:
                value = row
        return value

    def min_value(self, column):
        return self.find_value(column, min)
    def min_row(self, column):
        return self.find_row(column, min)
    def max_value(self, column):
        return self.find_value(column, max)
    def max_row(self, column):
        return self.find_row(column, max)

    def get_row(self, column, value):
        if not self.initialized:
            self.init()

        for row in self.rows:
            if row[column] == value:
                return row

        return None

    def get_rows(self, column, value):
        if not self.initialized:
            self.init()

        result = []
        for row in self.rows:
            if row[column] == value:
                result.append(row)

        return result

    def get_column(self, column):
        if not self.initialized:
            self.init()
        results = []
        for row in self.rows:
            results.append(row[column])
        return results

if __name__=='__main__':
    import doctest
    doctest.testmod()
