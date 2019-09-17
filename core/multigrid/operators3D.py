import numpy as np
from scipy import sparse
import itertools as iter

def set_laplacian(grid):
    """ 3D Laplacian with Neumann boundary conditions
    """
    nz, ny, nx = grid['shape']
#    nh = grid['nh']
    N = grid['N']
    G = grid['G']

     # get list of interior grid points
    idx = np.where(G > -1)
    # i the vector of 'x' coordinates
    # j the vector of 'y' coordinates
    k, j, i = idx[0], idx[1], idx[2]
    # the Laplacian matrix is N x N
    if N != len(i):
        raise ValueError('N should be equal to the number of domain grid points')

    # regular neighbours West, East, South and North
    listneighbours = [(-1, 0, 0), (0, -1, 0), (0, 0, -1), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    coefneighbours = [1]*6

    nmax = 7*N
    data = np.zeros((nmax,))
    rows = np.zeros((nmax,), dtype=int)
    cols = np.zeros((nmax,), dtype=int)

    count = 0
    for m in range(N):
        nbneighb = 0
        for who, coef in zip(listneighbours, coefneighbours):
            kk, jj, ii = who[0], who[1], who[2]
            inside = (
                ((i[m]+ii) >= 0) and
                ((i[m]+ii) <= nx-1) and
                ((j[m]+jj) >= 0) and
                ((j[m]+jj) <= ny-1) and
                ((k[m]+kk) >= 0) and
                ((k[m]+kk) <= nz-1))
            # set the coefficient only if the neighbour is inside the grid
            if inside:
                l = G[k[m]+kk, j[m]+jj, i[m]+ii]
                # and is interior
                if l > -1:
                    rows[count], cols[count], data[count] = m, l, coef
                    count += 1
                    nbneighb += coef

        # diagonal coef
        cff = -nbneighb
        rows[count], cols[count], data[count] = m, m, cff
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
    nzc, nyc, nxc = coarse['shape']
    nzf, nyf, nxf = fine['shape']
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nf*8
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    for idx in range(Nf):
        i, j, k = idx % nxf, (idx // nxf) % nyf,  idx // (nxf*nyf)
        idx2 = (i//2) + (j//2)*(nxc) + (k//2)*(nxc*nyc)
        di, cx = get_di_coef(i, nxf)
        dj, cy = get_di_coef(j, nyf)
        dk, cz = get_di_coef(k, nzf)
        data[l:l+8] = [a*b*c for a, b, c in iter.product(cz, cy, cx)]
        row[l:l+8] = [idx]*8
        col[l:l+8] = [idx2+c+b*nxc+a*nxc*nyc for a, b, c in iter.product(dk, dj, di)]
        l += 8
    I = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nf, Nc))
    return I.tocsr()

def set_coarsen(coarse, fine):
    nzc, nyc, nxc = coarse['shape']
    nzf, nyf, nxf = fine['shape']
    Nf = fine['N']
    Nc = coarse['N']
    nmax = Nc*8
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    l = 0
    print(coarse['shape'])
    for  idx in range(Nc):
        i, j, k = idx % nxc, (idx//nxc) % nyc, idx // (nxc*nyc)
        idx2 = (2*i) + (2*j)*(nxf) + (2*k)*(nxf*nyf)
        data[l:l+8] = [.25]*8
        row[l:l+8] = [idx]*8
        col[l:l+8] = [idx2, idx2+1, idx2+nxf, idx2+nxf+1,
                      idx2+nxf*nyf, idx2+nxf*nyf+1, idx2+nxf*nyf+nxf, idx2+nxf*nyf+nxf+1]
        l += 8
    R = sparse.coo_matrix((data[:l], (row[:l], col[:l])), shape=(Nc, Nf))
    return R.tocsr()
