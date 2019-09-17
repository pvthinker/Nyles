import numpy as np
import fortran_interp as fi


class Var(object): # to be renamed 'Scalar'
    def __init__(self, param):
        # nx = param.nx
        # ny = param.ny
        # nz = param.nz
        self.data = {}
        self.data['i'] = np.zeros((nx, nz, ny))
        self.data['j'] = np.zeros((ny, nx, nz))
        self.data['k'] = np.zeros((nz, ny, nx))
        self.activeview = 'k'

    def view(self, idx):
        """ return the 3D array
        copy it from the last activeview if necessary
        """
        
        if self.activeview == idx:
            field = self.data[idx]
        
        else:
            current = self.data[self.activeview]
            field = self.data[idx]
            if self.activeview+idx in ['ij', 'jk', 'ki']:
                rightswap(current, field)
                
            else:
                leftswap(current, field)

        self.activeview = idx
        return field
    
class Velocity(object): # to be renamed 'Vector'
    def __init__(self, param):
        self.u = Var(param)
        self.v = Var(param)
        self.w = Var(param)

    def view(self, idx):
        """ Return the three components of the velocity in the 'idx'
        outer index representation.
        """
        return [self.u.data[idx], self.v.data[idx], self.w.data[idx]] 

    def rswap(self, idx0, idx1):
        """ Copy the velocity from 'idx0' outer index representation
        to 'idx1' outer index representation, by swapping indices to
        the right, e.g.
        rswap('k', 'i') stands for 'kji' => 'ikj'

        idx0, idx1 = 'i', 'j' or 'k'
        """
        rightswap(self.u.data[idx0], self.u.data[idx1])
        rightswap(self.v.data[idx0], self.v.data[idx1])
        rightswap(self.w.data[idx0], self.w.data[idx1])
        
def rightswap(src, dest):
    dest[:, :, :] = np.transpose(src, (2, 0, 1))
#    dest[:, :, :] = np.moveaxis(src, 2, 0)

def leftswap(src, dest):
    dest[:, :, :] = np.transpose(src, (1, 2, 0))

def sxp(array2d):
    return np.roll(array2d, 1, axis=1) 

def syp(array2d):
    return np.roll(array2d, 1, axis=0) 

def sxm(array2d):
    return np.roll(array2d, -1, axis=1) 

def sym(array2d):
    return np.roll(array2d, -1, axis=0) 


def cov2contra(cova, contra):
        metric = 'sigma2D' # like ROMS/CROCO, h(x,y)
        metric = 'sigma1D' # h(x)
        #metric = 'cartesian' # dx, dy and dz are uniform, though not necessarily equal

        if metric == 'cartesian':
            u, v, w = cova.view('k')
            U, V, W = contra.view('k')
            U[:] = u*idx2
            V[:] = v*idy2
            W[:] = w*idz2

        elif metric == 'sigma1D':
            # cf Roullet et al, OM2017
            # slope = dz/dx, at cell center
            u, v, w = cova.view('j')
            U, V, W = contra.view('j')

            V[:] = v*idy2
            for j in range(ny):
                U[j][:, :] = u[j][:, :] - sxp(slope[j][:, :]*sym(w[j][:, :]))
                W[j][:, :] = gamma[j][:, :]*w[j][:, :] - syp(slope[j][:, :]*sxm(u[j][:, :]))
        else:
            pass

def div(cova, contra, divu):
    cov2contra(cova, contra)
    
    
def compute_vortexforce(vel, rhs, components='ijk'):
    """ Compute the vortex force in three steps:
    one step per component of the vorticity

    Makes the computation completely symmetrical with
    respect to components and index ordering

    After each step, duplicate the variables in
    arrays with reordered indices
    """

    if 'k' in components:
        # du += omega_z*v
        # dv -= omega_z*u
        # omega_z = dv/dx-du/dy
        u, v, w = vel.view('k')
        du, dv, dw = rhs.view('k')
        vortexforce(u, v, du, dv)
    vel.rswap('k', 'i')
    rhs.rswap('k', 'i')

    if 'j' in components:
        u, v, w = vel.view('j')
        du, dv, dw = rhs.view('j')
        vortexforce(v, w, dv, dw)
    vel.rswap('i', 'j')
    rhs.rswap('i', 'j')

    if 'i' in components:
        u, v, w = vel.view('i')
        du, dv, dw = rhs.view('i')
        vortexforce(w, u, dw, du)
    vel.rswap('j', 'k')
    rhs.rswap('j', 'k')

def addkengrav(ke, rho, rhs):
    """ add -grad(ke)-rho*grad(g*z)
    """
    pass

def rhs(state, t):
    vel = state[:3]
    p = state[3]
    rho = state[4]
    cov2contra(vel, contra)
    advect(rho, contra, drho)
    compute_vortexforce(vel, dvel, ke) # compute ke=u**2/2
    addkengrav(ke, rho, dvel) # add -grad(ke)-rho*grad(g*z)
    project(dvel, p)
    return dstate

def advect(trac, contra, dtrac):    
    for idx in 'ijk':
        phi = trac.view(idx)
        flux = contra.view(idx, component=idx)
        dphi = dtrac.view(idx)
        fo.upwind(phi, flux, dphi) # upwind transport in direction idx
    
def timescheme(state, t, dt):
    """ state is a list of scalars objects """
    if scheme == 'euler':
        dstate = rhs(state, t)
        for s, dsdt in zip(state, dstate):
            # s and dsdt are scalar objects
            # x and dxdt are 3D arrays
            x, dxdt = s.view('k'), dsdt.view('k')
            x += dxdt*dt

    elif scheme == 'RK3':
        pass
    
def vortexforce(uu, vv, duu, dvv, idx2=1.0, idy2=1.0):
    """ Add the omega_z component of the vortex force to the momentum

    du/dt += omega_z * V
    dv/dt -= omega_z * U

    where omega_z = ddx(v) - ddy(u)
    and (U,V) are the contravariant components of (u,v).
    If the grid is orthogonal, they read

    U = u/dx**2
    V = v/dy**2

    The routine assumes that the arrays are properly ordered,
    namely that
    - the outer index corresponds to the direction perpendicular
      to the (u,v) plane,
    - the middle index to the v direction,
    - and the inner index to the u direction.

    """
    nz, ny, nx = np.shape(uu)
    for k in range(nz): # outer loop
        u = uu[k]
        v = vv[k]
        du = duu[k]
        dv = dvv[k]
        if False:
            for j in range(1, ny-1):
                for i in range(1, nx-1):
                    # vorticity
                    omega = (v[j, i+1]-v[j, i])-(u[j+1, i]-u[j, i])
                    # v at u-point ! v is the contravariant component
                    vu = 0.25*(v[j, i]+v[j, i+1]+v[j-1, i]+v[j-1, i+1])
                    # u at v-point ! u is the contravariant component
                    uv = 0.25*(u[j, i]+u[j, i-1]+u[j+1, i]+u[j+1, i-1])

                    omega_j = omega # upwinded interpolation at u point along j
                    omega_i = omega # upwinded interpolation at v point along i

                    du[j, i] += vu * omega_j * idy2
                    dv[j, i] -= uv * omega_i * idx2
        elif False:
            omega = sxp(v)-v-syp(v)+u
            vu = 0.25*(v+sxp(v)+sym(v+sxp(v)))
            uv = 0.25*(u+syp(u)+sxm(u+syp(u)))
            du += vu*omega*idy2
            dv -= uv*omega*idx2
            
        else:
            fi.adv_upwind(u,v,du,dv,nh)

def assessperf(scalars):
    for dir in 'ijk':
        for p in scalars:
            p.view(dir)

nh = 3
nx=60+2*nh
ny=64+2*nh
nz=68+2*nh

param = {}

vel=Velocity(param)
rhs=Velocity(param)

idx2, idy2, idz2 = 1., 1., 1.
slope = Var(param)
gamma = Var(param)

scalars = []
for k in range(3):
    scalars.append(Var(param))

#a = %timeit -n 100 -o assessperf(scalars)
#print('best time : ', a.best)
