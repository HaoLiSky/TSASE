
import pyopencl as cl
import numpy as np

code=\
"""
__kernel void lj(__global float4* r, __global float* u, __global float4* f, unsigned int N)
{
    unsigned int i = get_global_id(0);
    unsigned int j = 0;
    float myu = 0.0;
    float fx = 0.0;
    float fy = 0.0;
    float fz = 0.0;
    for(j = 0; j < N; j++)
    {
        if(j!=i) // <== BAD MOJO
        {
            float dx = r[i].x - r[j].x;
            float dy = r[i].y - r[j].y;
            float dz = r[i].z - r[j].z;
            float d = sqrt(dx*dx + dy*dy + dz*dz);
            float invd = 1.0/d;
            float d3 = d * d * d;
            float d6 = d3 * d3;
            float d7 = d6 * d;
            float d12 = d6 * d6;
            float d13 = d12 * d;
            float fmag = (6.0/d7 - 12.0/d13);
            fx -= 4 * fmag * dx;
            fy -= 4 * fmag * dy;
            fz -= 4 * fmag * dz;
            myu += 0.5 * 4 * (1.0/d12 - 1.0/d6);
        }
    }

    u[i] = myu;
    f[i].x = fx;
    f[i].y = fy;
    f[i].z = fz;
}
"""

class ljocl:

    def __init__(self):
        self.atoms = None
        self.u = None
        self.f = None
        self.context = cl.create_some_context(False)
        self.queue = cl.CommandQueue(self.context)
        self.program = cl.Program(self.context, code).build()

    def get_potential_energy(self, atoms=None, force_consistent=False):
        if self.calculation_required(atoms, "energy"):
            self.atoms = atoms.copy()
            self.calculate()
        return self.u
        
    def get_forces(self, atoms):
        if self.calculation_required(atoms, "forces"):
            self.atoms = atoms.copy()
            self.calculate()
        return self.f.copy()
                        
    def get_stress(self, atoms):
        raise NotImplementedError
        
    def calculation_required(self, atoms, quantities):
        if atoms != self.atoms or self.atoms == None:
            return True
        if self.f == None or self.u == None or atoms == None:
            return True
        return False

    def set_atoms(self, atoms):
        pass

    def calculate(self):
        MF = cl.mem_flags
        r = self.atoms.get_positions()
        N = len(r)
        r = np.array(np.append(r,np.zeros(len(r)).reshape(len(r), 1), 1), dtype=np.float32)
        r = r.ravel().reshape(N, 4)
        u = np.zeros(N, dtype=np.float32)
        f = np.zeros((N,4), dtype=np.float32)
        rbuf = cl.Buffer(self.context, MF.READ_ONLY|MF.COPY_HOST_PTR, hostbuf = r)
        ubuf = cl.Buffer(self.context, MF.WRITE_ONLY, u.nbytes)
        fbuf = cl.Buffer(self.context, MF.WRITE_ONLY, f.nbytes)
        self.program.lj(self.queue, [N], None, rbuf, ubuf, fbuf, np.int32(N))
        cl.enqueue_read_buffer(self.queue, ubuf, u).wait()
        cl.enqueue_read_buffer(self.queue, fbuf, f).wait()
        rbuf.release()
        ubuf.release()
        fbuf.release()
        self.f = f
        self.u = u.sum()

                
        
        
