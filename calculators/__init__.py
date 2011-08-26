try:
    from al import al
except:
    pass
try:
    from morse import morse
except:
    pass
try:
    from lj import lj
except:
    pass
try:
    from ljocl import ljocl
except:
    pass
try:
    import lammps_xph
except:
    pass

def pt():
    return morse()
    

