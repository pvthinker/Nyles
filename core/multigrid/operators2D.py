import numpy as np
from scipy import sparse
import itertools as iter

def set_laplacian(grid):
    """ 2D Laplacian with Neumann boundary conditions
    """
    ny, nx = grid['shape']
#    nh = grid['nh']
    N = grid['N']
    G = grid['G']

     # get list of interior grid points
    idx = np.where(G > -1)
    # i the vector of 'x' coordinates
    # j the vector of 'y' coordinates
    j, i = idx[0], idx[1]
    # the Laplacian matrix is N x N
    if N != len(i):
        raise ValueError('N should be equal to the number of domain grid points')

    # regular neighbours West, East, South and North
    listneighbours = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    coefneighbours = [1]*4

    nmax = 5*N
    data = np.zeros((nmax,))
    rows = np.zeros((nmax,), dtype=int)
    cols = np.zeros((nmax,), dtype=int)

    count = 0
    for k in range(N):
        nbneighb = 0
        for who, coef in zip(listneighbours, coefneighbours):
            jj, ii = who[0], who[1]
            inside = (
                ((i[k]+ii) >= 0) and
                ((i[k]+ii) <= nx-1) and
                ((j[k]+jj) >= 0) and
                ((j[k]+jj) <= ny-1))
            # set the coefficient only if the neighbour is inside the grid
            if inside:
                l = G[j[k]+jj, i[k]+ii]
                # and is interior
                if l > -1:
                    rows[count], cols[count], data[count] = k, l, coef
                    count += 1
                    nbneighb += coef

        # diagonal coef
        cff = -nbneighb
        rows[count], cols[count], data[count] = k, k, cff
        count += 1

    A = sparse.coo_matrix(
        (data[:count], (rows[:count], cols[:count])), shape=(N, N))

    return A.tocsr()

def get_di_coef(i, nx):
    if i % 2 == 0:
        if i == 0:
            cx = [0, 1.]
            di = [0, 0]
        else:
            cx = [0.25, 0.75]
            di = [-1, 0]
    else:
        if i == nx-1:
            cx = [0, 1.]
            di = [0, 0]
        else:
            cx = [0.25, 0.75]
            di = [+1, 0]
    return di, cx
    
def set_interp(fine, coarse):
    nyc, nxc = coarse['shape']
    nyf, nxf = fine['shape']
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nf*4
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    for idx in range(Nf):
        i, j = idx % nxf, idx // nxf
        idx2 = (i//2) + (j//2)*(nxc)
        di, cx = get_di_coef(i, nxf)
        dj, cy = get_di_coef(j, nyf)
        data[l:l+4] = [a*b for a, b in iter.product(cy, cx)]
        row[l:l+4] = [idx]*4
        col[l:l+4] = [idx2+b+a*nxc for a, b in iter.product(dj, di)]
        l += 4
    I = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nf, Nc))
    return I.tocsr()

def set_coarsen(coarse, fine):
    nyc, nxc = coarse['shape']
    nyf, nxf = fine['shape']
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nc*4
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    for  idx in range(Nc):
        i, j = idx % nxc, idx//nxc
        idx2 = (2*i) + (2*j)*(nxf)
        data[l:l+4] = [.25]*4
        row[l:l+4] = [idx, idx, idx, idx]
        col[l:l+4] = [idx2, idx2+1, idx2+nxf, idx2+nxf+1]
        l += 4
    R = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nc, Nf))
    return R.tocsr()
