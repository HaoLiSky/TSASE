try:
    from lepspho import lepspho
except:
    pass
try:
    from al import al
except:
    pass
try:
    from morse import morse
    def pt():
        return morse()
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
    import lammps_ext
except:
    pass
try:
    from lisi import lisi
except:
    pass
try:
    from mo import mo
except:
    pass
try:
    from si import si
except:
    pass
try:
    from w import w
except:
    pass
try:
    from bopfox import bopfox
except:
    pass
try:
	from voter97 import voter97
except:
	pass
try:
	from ZDP_5Gauss import ZDP_5Gauss
except:
	pass

try:
    from socorro import Socorro
except:
    pass

try:
    from gauss3 import gauss3
except:
    pass


from push import push
from hyperplane import hyperplane_potential
