import ase.io
from con import read_con, write_con
from feff import read_feff, write_feff
from bopfox import read_bopfox, write_bopfox
from vasp import read_xdatcar
from socorro import read_socorro, write_socorro

def read(filename):
    try:
        return read_con(filename)
    except:
        pass
    try: 
        return read_bopfox(filename)
    except:
        pass
    try: 
        return read_xdatcar(filename)
    except:
        pass
    try:
        return read_socorro(filename)
    except:
        pass
    try:
        return ase.io.read(filename+"@:", format='vasp')
    except:
        pass
    try:
        return ase.io.read(filename+"@:")
    except:
        pass
    raise IOError, "Could not read file %s." % filename
