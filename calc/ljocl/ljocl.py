
import pyopencl as cl
import numpy as np

code=\
"""
#pragma OPENCL EXTENSION cl_amd_printf : enable
__kernel void lj(__global float4* r, __global float* u, __global float4* f, unsigned int N, float bx, float by, float bz)
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
            float dbx = dx/bx;
            float dby = dy/by;
            float dbz = dz/bz;
            float sx = dbx/fabs(dbx);
            float sy = dby/fabs(dby);
            float sz = dbz/fabs(dbz);
            dx = dx - bx * (int)(dbx + 0.5 * sx);
            dy = dy - by * (int)(dby + 0.5 * sy);
            dz = dz - bz * (int)(dbz + 0.5 * sz);
            float d = dx*dx + dy*dy + dz*dz;
            float invd = 1.0/d;
            float d2 = invd;
            float d6 = d2 * d2 * d2;
            float d12 = d6 * d6;
            float u = d12 - d6;
            float du = (12.0*d12 - 6.0*d6) * invd;
            fx += du * dx;
            fy += du * dy;
            fz += du * dz;
            myu += u;
        }
    }

    u[i] = 0.5 * 4 * myu;
    f[i].x = 4.0 * fx;
    f[i].y = 4.0 * fy;
    f[i].z = 4.0 * fz;
}
"""

class ljocl:

    def __init__(self):
        self.atoms = None
        self.u = None
        self.f = None
        self.ubuf = None
        self.fbuf = None
        self.rbuf = None
        self.lastSize = 0
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
        if self.lastSize != N:
            try:
                self.ubuf.release()
                self.fbuf.release()
                self.rbuf.release()
            except AttributeError:
                pass
            self.ubuf = cl.Buffer(self.context, MF.WRITE_ONLY, u.nbytes)
            self.fbuf = cl.Buffer(self.context, MF.WRITE_ONLY, f.nbytes)
            self.rbuf = cl.Buffer(self.context, MF.READ_ONLY, r.nbytes)
            self.lastSize = N
        cl.enqueue_write_buffer(self.queue, self.rbuf, r).wait()
        self.program.lj(self.queue, [N], None, self.rbuf, self.ubuf, self.fbuf, np.int32(N), 
                        np.float32(self.atoms.cell[0][0]),
                        np.float32(self.atoms.cell[1][1]),
                        np.float32(self.atoms.cell[2][2]))
        cl.enqueue_read_buffer(self.queue, self.ubuf, u)
        cl.enqueue_read_buffer(self.queue, self.fbuf, f)
        self.queue.finish()
        self.f = f
        self.u = u.sum()

                
        
        
