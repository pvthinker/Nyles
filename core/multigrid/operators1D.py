import numpy as np
from scipy import sparse

def set_laplacian(grid):
    """ 1D Laplacian with Neumann boundary conditions
    """
    nx = grid['shape'][0]
    nmax = 3*nx
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    for i in range(nx):
        idx = i
        left = (idx-1) % nx
        right = (idx+1) % nx
        row[l:l+3] = [idx, idx, idx]
        col[l:l+3] = [left, idx, right]
        cleft = 1*(i>0)
        cright = 1*(i<(nx-1))
        data[l:l+3] = [cleft, -(cleft+cright), cright]
        l += 3
    print('number of rows: %i' % nx)
    A = sparse.coo_matrix((data[:l], (row[:l], col[:])), shape=(nx, nx))
    return A.tocsr()

def set_interp(fine, coarse):
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nf*4
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    nx = Nf
    nxc = nx // 2
    for i in range(Nc):
        i1 = 2*i
        if i>0:
            cleft = 0.25
        else:
            cleft = 0.
        if i<(nxc-1):
            cright = 0.25
        else:
            cright = 0.
        data[l:l+4] = [cleft, 1-cleft, 1-cright, cright]
        row[l:l+4] = [i1, i1, (i1+1) % nx, (i1+1) % nx]
        col[l:l+4] = [(i-1) % nxc, i, i, (i+1) % nxc]
        l += 4
    I = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nf, Nc))
    return I.tocsr()

def set_coarsen(coarse, fine):
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
        data[l:l+2] = [.5, .5]
        row[l:l+2] = [i, i]
        col[l:l+2] = [i1, (i1+1) % nx]
        l += 2
    R = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nc, Nf))
    return R.tocsr()
