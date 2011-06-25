import numpy
import ase
import tsase
import time

head = "%5s %16s" % ("Natom", "ljocl speedup")
print head
print "=" * len(head)

ljocl = tsase.calc.ljocl()
lj = tsase.calc.lj()

for N in range(10, 400, 20):
    a = ase.Atoms(["H"]*N, positions=numpy.random.normal(0, 1, (N,3)))
    a.center(100)
    b = a.copy()
    a.set_calculator(ljocl)
    b.set_calculator(lj)

    Q = 10
    tcl = []
    t = []
    for i in range(Q):
        r = numpy.random.normal(0, 1, (N,3))
        a.set_positions(r.copy())
        b.set_positions(r.copy())
        t0 = time.time()
        fa = a.get_forces()
        tcl.append(time.time() - t0)
        t0 = time.time()
        fb = b.get_forces()
        t.append(time.time() - t0)
    print "%5d %16.12f" % (N, min(t)/min(tcl))

