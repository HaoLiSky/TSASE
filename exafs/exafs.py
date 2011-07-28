import sys
import os
import numpy
import ase
import math
import tsase
import tempfile
import subprocess
import shutil
import itertools
import atexit

__all__ = [ "load_chi_dat", "run_feff", "exafs" ]

def load_feff_dat(filename):
    k = []
    amplitude = []
    phase = []
    mean_free_path = []
    atoms = ase.Atoms()
    atoms.set_pbc((False,False,False))

    path_section = False
    atoms_section = False
    data_section = False
    f = open(filename)
    for line in f:
        line = line.strip()
        fields = line.split()
        if "genfmt" in line:
            path_section = True
            continue
        if fields[0] == "x" and fields[1] == "y" and fields[2] == "z":
            atoms_section = True
            path_section = False
            continue
        if fields[0] == "k" and fields[1] == "real[2*phc]":
            data_section = True
            atoms_section = False
            continue

        if path_section:
            if "---------------" in line:
                continue
            r_eff = float(fields[2])
            path_section = False

        if atoms_section:
            x = float(fields[0])
            y = float(fields[1])
            z = float(fields[2])
            pot = int(fields[3])
            atomic_number = int(fields[4])
            atoms.append(ase.Atom(symbol=atomic_number, position=(x,y,z),
                         tag=pot))

        if data_section:
            fields = [ float(f) for f in fields ]
            k.append(fields[0])
            amplitude.append(fields[2])
            phase.append(fields[3])
            mean_free_path.append(fields[5])

    k = numpy.array(k)
    amplitude = numpy.array(amplitude)
    phase = numpy.array(phase)
    mean_free_path = numpy.array(mean_free_path)
    return { 
             "atoms":atoms,
             "r_eff":r_eff,
             "k":k,
             "amplitude":amplitude,
             "phase":phase,
             "mean_free_path":mean_free_path,
           }

def load_files_dat(filename):
    begin = False
    filenames = []
    f = open(filename)
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        fields = line.split()
        if fields[0] == "file" and fields[1] == "sig2":
            begin = True
            continue

        if begin:
            filenames.append(fields[0])

    return filenames

def load_chi_dat(filename):
    f = open(filename)
    chi_section = False
    k = []
    chi = []
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue

        fields = line.split()
        if fields[0] == "k" and fields[1] == "chi" and fields[2] == "mag":
            chi_section = True
            continue

        if chi_section:
            k.append(float(fields[0]))
            chi.append(float(fields[1]))
    return numpy.array(k), numpy.array(chi)

def run_feff(atoms, absorber, tmp_dir=None):
    tmp_dir_path = tempfile.mkdtemp(prefix="tmp_feff_", dir=tmp_dir)
    tsase.io.write_feff(os.path.join(tmp_dir_path, "feff.inp"), atoms, absorber)
    p = subprocess.Popen(["feff"], cwd=tmp_dir_path, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    retval = p.wait()
    if retval != 0:
        print "Problem with feff calculation in %s" % tmp_dir_path
        return
    stdout, stderr = p.communicate()
    stderr = stderr.strip()
    if stderr == "hash error":
        atoms[absorber].set_position(atoms[absorber].get_position()+0.001)
        sys.stderr.write("%s\n"%stderr)
        return run_feff(atoms, absorber, tmp_dir)
    k, chi = load_chi_dat(os.path.join(tmp_dir_path, "chi.dat"))
    shutil.rmtree(tmp_dir_path)
    return k, chi

class DevNull:
    def write(self, string): pass
    def flush(self): pass
    def close(self): pass

def exafs(atoms, txt="-", tmp_dir=None, comm=None):
    from mpi4py import MPI

    if not comm:
        comm = MPI.COMM_WORLD

    size = comm.size
    rank = comm.rank

    if txt == "-" and rank == 0:
        outfile = sys.stdout
    elif txt == None and rank == 0:
        outfile = DevNull()
    elif rank == 0:
        outfile = open(txt, "a")
    else:
        outfile = DevNull()

    outfile.write("\nCalculating EXAFS Spectra\n")
    outfile.write("Number of Atoms: %i\n" % len(atoms))
    outfile.write("Number of Cores: %i\n" % size)
    outfile.flush()

    chi_total = {}
    atomic_symbols = set( [ a.symbol for a in atoms ] )
    for symbol in atomic_symbols:
        chi_total[symbol] = []

    k = None
    #absorbers = numpy.where((numpy.array(range(len(atoms)))%size)==rank)[0]
    
    for i in range(len(atoms)):
        if i%size != rank:
            continue
        k, chi = run_feff(atoms, i, tmp_dir)
        chi_total[atoms[i].symbol].append(chi)

    #in case more ranks than atoms
    k = comm.bcast(k)

    for symbol in atomic_symbols:
        if len(chi_total[symbol]) == 0:
            chi_total[symbol] = numpy.zeros(len(k))

        chi_total[symbol]  = numpy.sum(numpy.array(chi_total[symbol]), axis=0)
        chi_total[symbol]  = comm.allreduce(chi_total[symbol])
        natoms = len( [ a for a in atoms if a.symbol == symbol ] )
        chi_total[symbol] /= float(natoms)

    return k, chi_total

if __name__ == "__main__":
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.Get_rank()
    atoms = tsase.io.read_con(sys.argv[1])
    k, chi_total = exafs(atoms)

    if rank == 0:
        print k, chi_total
