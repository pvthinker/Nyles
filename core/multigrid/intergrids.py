"""

Provide the functions to coarsen and interpolate

"""

import numpy as np
import itertools as iter
import grids as gr
from mpi4py import MPI
from scipy import sparse
from timing import timing
import topology as topo
import gluetools


class Intergrids(object):
    def __init__(self, fine, coarse):
        self.fine = fine
        self.coarse = coarse
        # glueflag indicates if this intergrid needs glue/split
        self.glueflag = any(coarse.incr != fine.incr)
        if self.glueflag:
            self.glue = gluetools.Gluegrids(fine, coarse)
            self.dummy = self.glue.dummy
            self.Interpol = set_interpolation_matrix(fine, self.dummy)
            self.Restrict = set_restriction_matrix(fine, self.dummy)

        else:
            self.Interpol = set_interpolation_matrix(fine, coarse)
            self.Restrict = set_restriction_matrix(fine, coarse)

    @timing
    def fine2coarse(self, which='r'):
        assert which in "rb"

        y = self.dummy.x if self.glueflag else self.coarse.b
        y[:] = self.Restrict * self.fine.r

        if self.glueflag:
            # glue dummy.b onto coarse.b
            self.glue.glue_array("b")

        self.coarse.halofill('b')

        if which is 'b':
            self.coarse.r[:] = self.coarse.b

        self.coarse.x[:] = 0.

    @timing
    def coarse2fine(self):
        if self.glueflag:
            self.glue.split_array("x")
            x = self.dummy.x
        else:
            x = self.coarse.x

        self.fine.x += self.Interpol*x
        # halofill might be skept because
        # Interpol does the halo
        self.fine.halofill('x')


def get_di_coef(i, n, i0):
    """ 
    Provide the index/indices (di) and coef/coefs (cx)
    for the interpolation at point 'i' in the 'i' direction
    """
    if (i == 0) or (i == n-1):
        cx = [1.]
        di = [0]
    else:
        direc = -1 if (i-i0) % 2 == 0 else +1
        cx = [0.25, 0.75]
        di = [direc, 0]

    # if (i-i0) % 2 == 0:
    #     if (i == 0) or (i==n-1):
    #         cx = [1.]
    #         di = [0]
    #     else:
    #         cx = [0.25, 0.75]
    #         di = [-1, 0]
    # else:
    #     if (i==0) or (i == n-1):
    #         cx = [1.]
    #         di = [0]
    #     else:
    #         cx = [0.25, 0.75]
    #         di = [+1, 0]
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
    nh = 1  # coarse.nh
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
    K0, K1, J0, J1, I0, I1 = coarse.domainindices
    for k in range(k0, k1):

        if interpk:
            dk, ck = get_di_coef(k, nk, k0)
            k2 = K0+(k-k0)//2
        else:
            dk, ck = [0], [1.]
            k2 = k
        for j in range(j0, j1):
            if interpj:
                dj, cj = get_di_coef(j, nj, j0)
                j2 = J0+(j-j0)//2
            else:
                dj, cj = [0], [1.]
                j2 = j
            for i in range(i0, i1):
                if interpi:
                    di, ci = get_di_coef(i, ni, i0)
                    i2 = I0+(i-i0)//2
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

    k0, k1, j0, j1, i0, i1 = fine.domainindices
    K0, K1, J0, J1, I0, I1 = coarse.domainindices
    for k in range(K0, K1):
        if interpk:
            dk, ck = [0, 1], [0.5, 0.5]
            k2 = k0+(k-K0)*2
        else:
            dk, ck = [0], [1.]
            k2 = k
        for j in range(J0, J1):
            if interpj:
                dj, cj = [0, 1], [0.5, 0.5]
                j2 = j0+(j-J0)*2
            else:
                dj, cj = [0], [1.]
                j2 = j
            for i in range(I0, I1):
                if interpi:
                    di, ci = [0, 1], [0.5, 0.5]
                    i2 = i0+(i-I0)*2
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

    procs = [4, 4, 1]
    topo.topology = 'perio_y'
    myrank = 0
    nh = 1
    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs, extension=26)

    print('-'*80)
    nx, ny, nz = 8, 8, 1
#    nx, ny, nz = 4, 4, 4

    param = {'nx': nx, 'ny': ny, 'nz': nz, 'nh': nh, 'procs': procs,
             "nglue": 32, "ncellscoarsest": 16,
             "omega": 0.8, "ndeepest": 32}

    allgrids = gr.define_grids(param)

    nlevs = len(allgrids)

    subdomains = subdom.set_subdomains(allgrids)
    subdom.attach_subdomain_to_grids(allgrids, subdomains, myrank)
    # subdom.print_subdomains(subdomains)

    grids = []
    for lev in range(nlevs):
        grd = allgrids[lev]
        subdomain = subdomains[grd["subdomain"]]
        subdomain['myrank'] = myrank
        grd['neighbours'] = subdomain['allneighbours'][myrank]
        grd['extension'] = 26
        grids += [gr.Grid(grd, param)]

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
    print(coarse.infos["glueflag"])
    print(coarse.infos["subd"])
