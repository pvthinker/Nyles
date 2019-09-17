import numpy as np
from scipy import sparse
#from scipy.sparse import linalg as LA
from numpy import linalg as LA


class Multigrid(object):
    def __init__(self, param):
        
        self.grid = self.guess_hierarchy()
        self.npre = 4
        self.npost = 4
        self.ndeepest = 10
        self.verbose = True
        self.tol = 1e-12
        self.maxite = 4
        self.omega = .85 # for the smoother
        self.nlevs = len(self.grid)
        self.setup_grids()

    def guess_hierarchy(self):
        nx, ny, nz = 64, 96, 96
        nlevs = 0
        grids = []
        done = False
        while not(done):
            shape = (nz, ny, nx)
            N = np.asarray(shape).prod()
            g = {'shape': shape, 'N': N}
            grids.append(g)
            if (nx % 2 == 0) and (nx > 3):
                nx = nx//2
                ny = max(1, ny//2)
                nz = max(1, nz//2)
            else:
                if nx <= 3:
                    done = True
                else:
                    raise ValueError('pb with nx')
            nlevs += 1
        print('Levels are :')
        for g in grids:
            print(g['shape'])
        return grids
        
    def setup_grids(self):
        """ main initialisation part of the multigrid """
        # define the grids hierarchy
        self.x = []
        self.b = []
        self.r = []

        import operators3D as op
        # populate them
        for lev in range(self.nlevs):
            if lev < self.nlevs-1:
                N = self.grid[lev]['N']
                Nc = self.grid[lev+1]['N']
                # I = define_interp(self.grid[lev], self.grid[lev+1])
                I = op.set_interp(self.grid[lev], self.grid[lev+1])
                # data = np.ones((Nc,))
                # I = sparse.spdiags(data, 0, N, Nc)
                # I = I.tocsr()
                self.grid[lev]['I'] = I
                # self.grid[lev]['R'] = I.T *.5
                self.grid[lev]['R'] = op.set_coarsen(self.grid[lev+1], self.grid[lev])
#                self.grid[lev]['I'].sort_indices()
#                self.grid[lev]['R'].sort_indices()

            if lev == 0:
                N = self.grid[lev]['N']
                shape = self.grid[lev]['shape']
#                A = define_smoother(self.grid[lev])
                G = np.reshape(np.arange(N, dtype=int), shape)
                self.grid[lev]['G'] = G
                
                A = op.set_laplacian(self.grid[lev])

            else:
                I = self.grid[lev-1]['I']
                R = self.grid[lev-1]['R']
                Af = self.grid[lev-1]['A']
                A = R * Af * I

            D = A.diagonal()
            self.grid[lev]['S'] = A-sparse.diags(D)
            if all(D == D[0]):
                self.grid[lev]['idiag'] = -1./D[0]
            else:
                self.grid[lev]['idiag'] = -1./D
            
            self.grid[lev]['A'] = A
#            self.grid[lev]['A'].sort_indices()
#            self.grid[lev]['S'].sort_indices()
                
            N = self.grid[lev]['N']
            
            self.x += [np.zeros((N,))]
            self.b += [np.zeros((N,))]
            self.r += [np.zeros((N,))]

    def solve(self, x, b, cycle='V'):
        """ driver to solve A*x=b using the Fcycle
        x is the first guess, b the right hand side"""
        
        self.x[0][:] = x
        self.b[0][:] = b
        
        g = self.grid[0]
        self.residual(0)
        normeb = self.norm(b)
        if self.verbose:
            print('||b|| = ', normeb)
            
        if normeb > 0:
            res0 = self.norm(self.r[0])/normeb
        else:
            return 0, 0.
        res = res0
        nite = 0
        nite_diverge = 0
        # improve the solution until one of this condition is wrong
        while (nite < self.maxite) and (res0 > self.tol):
            if cycle == 'V':
                self.Vcycle(0)
            elif cycle == 'F':
                self.Fcycle()
            else:
                raise ValueError('use cycle V or F')
            
            self.residual(0)

            res = self.norm(self.r[0])/normeb
            conv = res0/res

            res0 = res
            nite += 1
            if self.verbose:
                print(' ite = {} / res = {:.2e} / conv = {:8.4f}'.format(nite, res, conv))

            if (conv < 1):
                nite_diverge += 1

            if (nite_diverge > 4):
                raise ValueError('solver is not converging')

        x[:] = self.x[0]
        return nite, res

    def Fcycle(self):
        for lev in range(self.nlevs-1):
            if lev == 0:
                r = self.r[0]
            else:
                r = self.b[lev]
            self.fine2coarse(lev, r)

        self.solveexactly()

        for lev in range(self.nlevs-2, -1, -1):
            self.coarse2fine(lev)
            self.Vcycle(lev)

    def twolev(self, lev):
        self.smooth(lev, self.npre)
        self.residual(lev)
        self.fine2coarse(lev, self.r[lev])
        self.smooth(lev+1, 500)
        self.coarse2fine(lev)
        self.smooth(lev, self.npost)
        
    def Vcycle(self, lev0):
        for lev in range(lev0, self.nlevs-1):
            self.smooth(lev, self.npre)
            self.residual(lev)
            self.fine2coarse(lev, self.r[lev])

        self.solveexactly()

        for lev in range(self.nlevs-2, lev0-1, -1):
            self.coarse2fine(lev)
            self.smooth(lev, self.npost)
    
    def smooth(self, lev, nite):
        for k in range(nite):
            self.x[lev][:] = self.x[lev]*(1-self.omega) + (self.grid[lev]['S']*self.x[lev]-self.b[lev])*(self.omega*self.grid[lev]['idiag'])

    def resnorm(self, lev):
        self.residual(lev)
        return self.norm(self.r[lev])
    
    def residual(self, lev):
        self.r[lev][:] = self.b[lev]-self.grid[lev]['A']*self.x[lev]

    def norm(self, phi):
        return LA.norm(phi)
        
    def solveexactly(self):
        lev = -1 # last level = coarsest grid
        self.smooth(lev, self.ndeepest)
        #self.x[lev][:] = LA.solve(self.grid[lev]['A'], self.b[lev])
    
    def coarse2fine(self, lev):
        # interpolate grid(lev+1)%p to grid(lev)%r and add it to grid(lev)%p
        self.x[lev] += self.grid[lev]['I']*self.x[lev+1]
        #print(lev, 'mean(x)=', np.mean(self.x[lev]))

    def fine2coarse(self, lev, r):
        # coarsen grid(lev)%r to grid(lev+1)%b and set grid(lev+1)%p=0
        self.b[lev+1][:] = self.grid[lev]['R']*r
        self.x[lev+1][:] = 0.

def define_smoother(grid):
    nx = grid['shape'][0]
    nmax = 3*nx
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    for i in range(nx):
        idx = i
        row[l:l+3] = [idx, idx, idx]
        col[l:l+3] = [(idx-1) % nx, idx, (idx+1) % nx]
        data[l:l+3] = [1, -2, 1]
        l += 3
    print('number of rows: %i' % nx)
    A = sparse.coo_matrix((data[:l], (row[:l], col[:])), shape=(nx, nx))
    return A.tocsr()
       
def laplacian2D(grid):
    ny, nx = grid['shape']
    N = nx*ny
    nmax = 5*N
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    for j in range(ny):
        for i in range(nx):
            idx = i+j*ny
            im = (i-1) % nx + j*nx
            ip = (i+1) % nx + j*nx
            jm = i + ((j-1) % ny)*nx
            jp = i + ((j+1) % ny)*nx
            row[l:l+5] = [idx, idx, idx, idx, idx]
            col[l:l+5] = [jm, im, idx, ip, jp]
            data[l:l+5] = [1, 1, -4, 1, 1]
            l += 5
    print('number of rows: %i' % N)
    A = sparse.coo_matrix((data[:l], (row[:l], col[:])), shape=(N, N))
    return A.tocsr()
       
def define_smoother3(grid):
    nh = grid.nh
    nx = grid.nx
    ny = grid.ny
    nz = grid.nz
    cx = grid.cx
    cy = grid.cy
    cz = grid.cz

    lx = 1
    ly = (nx+2*nh)
    lz = (ny+2*nh)*ly
    l = 0
    nmax = 7*nx*ny*nz
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                idx = (i+nh)+((j+nh)+(k+nh)*lz)*ly
                diag = 0
                if i>0:
                    row[l], col[l], data[l] = idx, idx-lx, cx
                    diag += cx 
                    l += 1
                if i<nx-1:
                    row[l], col[l], data[l] = idx, idx+lx, cx
                    diag += cx 
                    l += 1
                if j>0:
                    row[l], col[l], data[l] = idx, idx-ly, cy
                    diag += cy 
                    l += 1
                if j<ny-1:
                    row[l], col[l], data[l] = idx, idx+ly, cy
                    diag += cy 
                    l += 1
                if k>0:
                    row[l], col[l], data[l] = idx, idx-lz, cz
                    diag += cz 
                    l += 1
                if k<nz-1:
                    row[l], col[l], data[l] = idx, idx+lz, cz
                    diag += cz
                    l += 1
                if diag>0:
                    row[l], col[l], data[l] = idx, idx, -diag
                    l += 1
    print('number of rows: %i' % l)
    A = sparse.coo_matrix((data[:l], (row[:l], col[:])), shape=(l, l))
    return A.tocsr()

def define_interp4(fine, coarse):
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nf*4
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    nx = Nf
    for i in range(Nc):
        i1 = 2*i
        data[l:l+4] = [0.25, 0.75, 0.75, 0.25]
        row[l:l+4] = [(i1-1) % nx, i1, (i1+1) % nx, (i1+2) % nx]
        col[l:l+4] = [i, i, i, i]
        l += 4
    I = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nf, Nc))
    return I.tocsr()
        
def define_interp2(fine, coarse):
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nf*2
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    nx = Nf
    for i in range(Nc):
        i1 = 2*i
        data[l:l+2] = [1., 1.]
        row[l:l+2] = [(i1-1) % nx, i1]
        col[l:l+2] = [i, i]
        l += 2
    I = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nf, Nc))
    return I.tocsr()
        
def define_interp(fine, coarse):
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nf*3
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    nx = Nf
    for i in range(Nc):
        i1 = 2*i
        data[l:l+3] = [0.5, 1., 0.5]
        row[l:l+3] = [(i1-1) % nx, i1, (i1+1) % nx]
        col[l:l+3] = [i, i, i]
        l += 3
    I = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nf, Nc))
    return I.tocsr()
        
def define_interp2D(fine, coarse):
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nf*3
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    nx = Nf
    for i in range(Nc):
        i1 = 2*i
        data[l:l+3] = [0.5, 1., 0.5]
        row[l:l+3] = [(i1-1) % nx, i1, (i1+1) % nx]
        col[l:l+3] = [i, i, i]
        l += 3
    I = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nf, Nc))
    return I.tocsr()
        
     
# def define_interpolation(fine, coarse):
#     for k in range(nz):
#         for j in range(ny):
#             for i in range(nx):
#                 idx = (i+nh)+((j+nh)+(k+nh)*lz)*ly

#                 row[l:l+8] = [

param = {}
nx = 64
m = Multigrid(param)
ndims = len(m.grid[0]['shape'])
if ndims == 1:
    b = np.zeros((nx,))
    i0 = nx//3
    b[i0] = 1
    b[i0+nx//2] = -1
    x = b*0
elif ndims == 2:
    b = np.zeros(m.grid[0]['shape'])
    b[10, 10] = 1.
    b[-10, -10] = -1.
    b = b.ravel()
    x = b*0.
elif ndims == 3:
    b = np.zeros(m.grid[0]['shape'])
    b[10, 10, 0] = 1.
    b[-10, -10, -1] = -1.
    b = b.ravel()
    x = b*0.
    
#m.solve(x,b)
m.b[0][:] = b
#m.twolev(0)
#x = m.x[0]

from matplotlib import pyplot as plt
plt.ion()

if False:
    for k in range(5):
        plt.clf()
        m.residual(0)
        print('r0=',m.resnorm(0))
        plt.plot(m.r[0],label='r0')
        m.fine2coarse(0, m.r[0])

        plt.plot(m.b[1],label='b1')
        m.smooth(1, 500)
        plt.plot(m.x[1],label='x1')
        print('r1=',m.resnorm(1))
        m.coarse2fine(0)
        m.smooth(0, 5)
        plt.plot(m.x[0],label='x0')
        print('r0=',m.resnorm(0))
        plt.plot(m.r[0],alpha=0.5,label='r0')
        plt.legend()

if False:
    plt.clf()
    m.x[0][:] = 0
    m.b[0][:] = b
    for k in range(10):
        m.residual(0)
        print('r0=',m.resnorm(0))
        m.Fcycle()
        plt.plot(m.x[0],label='%i' % k)
    plt.legend()
