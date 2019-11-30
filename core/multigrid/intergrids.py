"""

Provide the functions to coarsen and interpolate

"""

import numpy as np
import itertools as iter
from scipy import sparse
from timing import timing

class Intergrids(object):
    def __init__(self, fine, coarse):
        self.fine = fine
        self.coarse = coarse
        self.Interpol = set_interpolation_matrix(fine, coarse)
        self.Restrict = set_restriction_matrix(fine, coarse)
        # todo: implement glue and split
        self.gluing = False
        if self.gluing:
            # define self.dummy and gluing process
            raise ValueError('not yet implemented')

    @timing
    def fine2coarse(self, which='r'):
        assert which in ['r', 'b']

        y = self.dummy if self.gluing else self.coarse.b
        y[:] = self.Restrict * self.fine.r

        if self.gluing:
            self.glue(self.dummy, self.coarse.b)

        self.coarse.halofill('b')

        if which is 'b':
            self.coarse.r[:] = self.coarse.b

        self.coarse.x[:] = 0.

    @timing
    def coarse2fine(self):
        if self.gluing:
            self.split(self.coarse.x, self.dummy)
            x = self.dummy
        else:
            x = self.coarse.x

        self.fine.x += self.Interpol*x
        # halofill might be skept because
        # Interpol does the halo
        # self.fine.halofill('x')

    def glue(self, x, y):
        """ glue x's into y """
        pass

    def split(self, x, y):
        """ split x into y """
        pass


def get_di_coef(i, n, i0):
    """ 
    Provide the index/indices (di) and coef/coefs (cx)
    for the interpolation at point 'i' in the 'i' direction
    """
    if (i-i0) % 2 == 0:
        if (i == 0):
            cx = [1.]
            di = [0]
        else:
            cx = [0.25, 0.75]
            di = [-1, 0]
    else:
        if (i == n-1):
            cx = [1.]
            di = [0]
        else:
            cx = [0.25, 0.75]
            di = [+1, 0]
    return di, cx


def set_interpolation_matrix(fine, coarse):
    """

    Define Interpolation matrix

    linear interpolation with (0.25, 0.75) weights

    near boundary use 1. weight

    function is general enough to handle interpolation
    in 1D, 2D or 3D. Direction to interpolate are
    recovered by comparing shapes of 'fine' and 'coarse'

    Interpolation should be performed in the halo (!)
    not to by-pass halo fill but to define the Laplacian
    on the coarse grid (see mg.py for more explanations)

    """
    nh = coarse.nh
    nzc, nyc, nxc = coarse.shape
    nzf, nyf, nxf = fine.shape
    # determine which direction has to be interpolated
    interpk = (nzf//2 == nzc)
    interpj = (nyf//2 == nyc)
    interpi = (nxf//2 == nxc)

    Nf = fine.N
    Nc = coarse.N
    nmax = Nf*8
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    counter = 0

    nk, nj, ni = fine.size

    Gc = np.arange(coarse.N).reshape(coarse.size)
    Gf = np.arange(fine.N).reshape(fine.size)

    k0, k1, j0, j1, i0, i1 = fine.domainindices
    for k in range(nk):

        if interpk:
            dk, ck = get_di_coef(k, nk, k0)
            k2 = k0+(k-k0)//2
        else:
            dk, ck = [0], [1.]
            k2 = k
        for j in range(nj):
            if interpj:
                dj, cj = get_di_coef(j, nj, j0)
                j2 = j0+(j-j0)//2
            else:
                dj, cj = [0], [1.]
                j2 = j
            for i in range(ni):
                if interpi:
                    di, ci = get_di_coef(i, ni, i0)
                    i2 = i0+(i-i0)//2
                else:
                    di, ci = [0], [1.]
                    i2 = i

                Js = [Gc[k2+a, j2+b, i2+c]
                      for a, b, c in iter.product(dk, dj, di)]
                nblock = len(Js)
                Is = [Gf[k, j, i]]*nblock
                coefs = [a*b*c for a, b, c in iter.product(ck, cj, ci)]
                data[counter:counter+nblock] = coefs
                row[counter:counter+nblock] = Is
                col[counter:counter+nblock] = Js
                counter += nblock
    I = sparse.coo_matrix((data[:counter], (row[:counter], col[:counter])),
                          shape=(Nf, Nc))
    return I.tocsr()


def set_restriction_matrix(fine, coarse):
    """
    Define the restriction matrix

    restriction is defined by a basic two points averaging
    
    restriction can be 1D, 2D or 3D. 
    """
    nh = coarse.nh
    nzc, nyc, nxc = coarse.shape
    nzf, nyf, nxf = fine.shape
    interpk = (nzf//2 == nzc)
    interpj = (nyf//2 == nyc)
    interpi = (nxf//2 == nxc)

    Nf = fine.N
    Nc = coarse.N
    nmax = Nc*8
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    counter = 0

    Gc = np.arange(coarse.N).reshape(coarse.size)
    Gf = np.arange(fine.N).reshape(fine.size)

    k0, k1, j0, j1, i0, i1 = coarse.domainindices
    for k in range(k0, k1):
        if interpk:
            dk, ck = [0, 1], [0.5, 0.5]
            k2 = k0+(k-k0)*2
        else:
            dk, ck = [0], [1.]
            k2 = k
        for j in range(j0, j1):
            if interpj:
                dj, cj = [0, 1], [0.5, 0.5]
                j2 = j0+(j-j0)*2
            else:
                dj, cj = [0], [1.]
                j2 = j
            for i in range(i0, i1):
                if interpi:
                    di, ci = [0, 1], [0.5, 0.5]
                    i2 = i0+(i-i0)*2
                else:
                    di, ci = [0], [1.]
                    i2 = i

                Js = [Gf[k2+a, j2+b, i2+c]
                      for a, b, c in iter.product(dk, dj, di)]
                nblock = len(Js)
                Is = [Gc[k, j, i]]*nblock
                coefs = [a*b*c for a, b, c in iter.product(ck, cj, ci)]
                data[counter:counter+nblock] = coefs
                row[counter:counter+nblock] = Is
                col[counter:counter+nblock] = Js
                counter += nblock
    R = sparse.coo_matrix((data[:counter], (row[:counter], col[:counter])),
                          shape=(Nc, Nf))
    return R.tocsr()


# ----------------------------------------------------------------------
if __name__ == '__main__':

    import grids as gr
    import topology as topo
    import subdomains as subdom

    procs = [1, 1, 1]
    topo.topology = 'perio_y'
    myrank = 0
    nh = 1
    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs, extension=26)

    print('-'*80)
    nx, ny, nz = 8, 8, 1
#    nx, ny, nz = 4, 4, 4

    param = {'nx': nx, 'ny': ny, 'nz': nz, 'nh': nh, 'procs': procs}

    allgrids = gr.define_grids(param)

    nlevs = len(allgrids)

    subdomains = subdom.set_subdomains(allgrids)
    # subdom.print_subdomains(subdomains)

    grids = []
    for lev in range(nlevs):
        grd = allgrids[lev]
        subdomain = subdom.isolate(subdomains, lev, myrank)
        subdomain['myrank'] = myrank
        grd['neighbours'] = subdomain['neighbours']
        grd['extension'] = 26
        grids += [grd.Grid(grd)]

    lev = 0
    fine = grids[lev]
    coarse = grids[lev+1]

    I = set_interpolation_matrix(fine, coarse)
    R = set_restriction_matrix(fine, coarse)

    print('****** interpolation')
    coarse.x[:] = np.arange(coarse.N)
    fine.x[:] = I*coarse.x
    # fine.halofill('x')
    print(coarse.toarray('x'))
    print(fine.toarray('x'))

    print('****** restriction')
    fine.x[:] = np.arange(fine.N)
    coarse.x[:] = R*fine.x
    print(fine.toarray('x'))
    print(coarse.toarray('x'))

    print('****** Intergrids')
    intg = Intergrids(fine, coarse)
    fine.r[:] = np.arange(fine.N)
    intg.fine2coarse(which='b')
    print(coarse.toarray('b'))

    fine.x[:] = 0
    coarse.x[:] = np.arange(coarse.N)
    intg.coarse2fine()
    print(fine.toarray('x'))

    print(coarse.nh)
