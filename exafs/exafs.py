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

#def calculate_chi(feff):
#    chi = numpy.zeros(len(feff[0]["k"]))
#    for shell in feff:
#        for i in range(len(shell["k"])):
#            k = shell["k"][i]
#            f = shell["amplitude"][i]
#            d = shell["phase"][i]
#            r = shell["r_eff"]
#            l = shell["mean_free_path"][i]
#
#            if k==0.0:
#                continue
#            chi_partial  = f*(len(shell["atoms"])-1)/(k*r*r)
#            chi_partial *= math.exp(-2*r/l)*math.sin(2*k*r+d)
#            chi[i] += chi_partial
#    return feff[0]["k"],chi

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

def run_feff(atoms, absorber, tmp_dir="."):
    tmp_dir_path = tempfile.mkdtemp(prefix="tmp_feff_", dir=tmp_dir)
    tsase.io.write_feff(os.path.join(tmp_dir_path, "feff.inp"), atoms, absorber)
    p = subprocess.Popen(["feff"], cwd=tmp_dir_path, stdout=subprocess.PIPE)
    retval = p.wait()
    if retval != 0:
        print "Problem with feff calculation in %s" % tmp_dir_path
        return
    k, chi = load_chi_dat(os.path.join(tmp_dir_path, "chi.dat"))
    shutil.rmtree(tmp_dir_path)
    return k, chi

def exafs(atoms, tmp_dir="."):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    chi_total = []

    for i in range(len(atoms)):
        if i%size != rank:
            continue
        k, chi = run_feff(atoms, i, tmp_dir)
        chi_total.append(chi)

    chi_total = numpy.average(numpy.array(chi_total), axis=0)
    chi_total = comm.gather(chi_total)

    if rank == 0:
        chi_total = numpy.array(chi_total)
        chi_total = numpy.average(chi_total, axis=0)
        return k,chi_total
    return (None, None)

if __name__ == "__main__":
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.Get_rank()
    atoms = tsase.io.read_con(sys.argv[1])
    k, chi_total = exafs(atoms)
    if rank == 0:
        import pylab
        pylab.plot(k,chi_total)
        pylab.show()
