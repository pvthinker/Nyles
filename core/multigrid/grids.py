"""

  grids.py provides functions to define the multigrid hierarchy

  the dimensions of each grid, the corresponding domain decomposition,
  the operations to be done to go from one level to another

  these operations are any combination of:

  coarsening in 'i', 'j' and/or 'k',
      gluing in 'i', 'j' and/or 'k'

  coarsening halves the number of grid cells in the corresponding direction
  gluing halves the number of subdomains in the corresponding direction

  on the finest grid (level=0), each core computes its own subdomain
  on the coarsest grid, all cores compute the whole domain

  if you use this module as a main, it is suggested to redirect the
  output in a text file as the output is quite lengthy

  % python grids.py > out.txt

"""
import numpy as np
import halo
import fortran_operators as fortran
from scipy import sparse
from timing import timing
import mpitools
from mpi4py import MPI
import topology as topo


class Grid(object):
    """

    Provide functions working on a particular grid

      - norm
      - smooth
      - residual
      - solveexactly

    """

    def __init__(self, grid, param):
        for key in ['shape', 'incr', 'procs']:
            setattr(self, key, grid[key])

        for key in ['omega', 'ndeepest', 'myrank']:
            setattr(self, key, param[key])

        nh = 1
        procs0 = [i*p for i, p in zip(self.incr, self.procs)]
        print("---------> procs0=", procs0, self.myrank)
        loc = topo.rank2loc(self.myrank, procs0)
        neighbours = topo.get_neighbours(loc, procs0, incr=self.incr)
        size, domainindices = topo.get_variable_shape(
            self.shape, neighbours, nh)
        N = np.prod(size)

        self.domainindices = domainindices
        self.size = size
        self.N = N
        self.nh = nh
        self.neighbours = neighbours

        self.x = np.zeros((self.N,))
        self.b = np.zeros((self.N,))
        self.r = np.zeros((self.N,))
        self.halo = halo.Halo({"nh": nh, "size": size, "neighbours": neighbours,
                               "domainindices": domainindices, "shape": self.shape})

    def toarray(self, x):
        """
        Return the 3D array view of 'x'
        """
        assert x in ['x', 'b', 'r']
        return np.reshape(getattr(self, x), self.size)

    @timing
    def halofill(self, x):
        """
        Fill x's halo, where 'x' is in ['xbr']
        """
        self.halo.fill(self.toarray(x))

    def set_ADS(self, A):
        self.A = A

        D = A.diagonal()
        self.S = A-sparse.diags(D)

        idx = np.where(D != 0)
        self.idiag = np.zeros_like(D)
        self.idiag[idx] = -1./D[idx]

    @timing
    def norm(self, which='r'):
        """
        Compute the norm of 'which'
        excluding the halo
        """
        y = self.toarray(which)
        k0, k1, j0, j1, i0, i1 = self.domainindices
        localsum = fortran.norm(y, k0, k1, j0, j1, i0, i1)
        # Note: the global sum is done on  all cores,
        # meaning it works only on the finest partition,
        # viz. one subdomain per core.
        # BUT since the norm is needed only on the finest level,
        # we don't need to generalize this sum to
        # other partitions
        y2 = mpitools.global_sum(localsum)
        return np.sqrt(y2)

    @timing
    def residual(g):
        # g is self
        g.r[:] = g.b - g.A * g.x
        g.halofill('r')

    @timing
    def smooth(g, nite):
        # g is self
        for k in range(nite):
            g.x[:] = g.x*(1.-g.omega) + (g.S*g.x-g.b)*(g.omega*g.idiag)
            g.halofill('x')

    def solveexactly(g):
        # g is self
        g.smooth(g.ndeepest)


def define_grids(param):
    """
    define the grid hierarchy along with the coarsening/gluing operations to do

    gluing is triggered when the size in that direction is <= 'nglue'

    a grid is declared the coarsest when all cores compute the whole domain
    and that the number of cells of the global domain is <= 'ncellscoarsest'
    """
    procs = param['procs'].copy()
    nx, ny, nz = param['nx'], param['ny'], param['nz']
    nh = 1  # param['nh']

    nglue = param['nglue']
    ncellscoarsest = param['ncellscoarsest']

    shape = [nz, ny, nx]
    lev = 0
    done = False

    grid = {'nh': nh, 'shape': shape.copy(), 'procs': procs.copy(),
            'coarsen': {}}
    mygrids = [grid]

    while not(done):
        r = {}
        done = True
        nglued = ''
        for k, direc in enumerate('kji'):
            if ((shape[k] % 2) == 0):
                r['r'+direc] = 1
                shape[k] //= 2
                done = False
            elif (shape[k] > 3):
                raise ValueError(
                    'size should be 3*2**p or 2**p in direction %s' % direc)
            else:
                pass
            if (shape[k] <= nglue) and (procs[k] > 1):
                shape[k] *= 2
                procs[k] //= 2
                r['g'+direc] = 1
                nglued += direc

        if len(nglued) == 1:
            # if there is one gluing, let's do another one
            nprocs = np.prod(procs)
            # print('only one glued, nprocs=%i' % nprocs)
            if nprocs == 1:
                # there is nothing we can do
                pass
            elif nprocs >= 2:
                direc = nglued
                k = 'kji'.index(direc)
                if procs[k] > 1:
                    shape[k] *= 2
                    procs[k] //= 2
                    r['g'+direc] = 2
                    # print('glue two instead of one')
                else:
                    # wait another turn
                    shape[k] //= 2
                    procs[k] *= 2
                    r.pop('g'+direc)
                    # print('postpone gluing to next level')

        grid = {'nh': nh, 'shape': shape.copy(), 'procs': procs.copy(),
                'coarsen': r.copy()}
        mygrids += [grid]
        lev += 1
        if np.prod(shape) <= ncellscoarsest:
            # assert np.prod(procs) == 1, 'weird problem with the grid hierachy'
            done = True

    return mygrids


def define_grids_v2(param):
    """
    define the grid hierarchy along with the coarsening/gluing operations to do

    gluing is triggered when the size in that direction is <= 'nglue'

    a grid is declared the coarsest when all cores compute the whole domain
    and that the number of cells of the global domain is <= 'ncellscoarsest'
    """
    procs = param['procs']
    incr = np.array([1, 1, 1], dtype=int)
    nx, ny, nz = param['nx'], param['ny'], param['nz']
    shape = [nz, ny, nx]

    nh = 1

    nglue = param['nglue']
    ncellscoarsest = param['ncellscoarsest']

    lev = 0
    done = False

    grid = {'shape': shape.copy(), 'procs': procs.copy(),
            'incr': incr.copy(), 'coarsen': {}}
    mygrids = [grid]

    while not(done):
        r = {}
        done = True
        nglued = ''
        for k, direc in enumerate('kji'):
            if ((shape[k] % 2) == 0):
                r['r'+direc] = 1
                shape[k] //= 2
                done = False
            elif (shape[k] > 3):
                raise ValueError(
                    'size should be 3*2**p or 2**p in direction %s' % direc)
            else:
                pass
            if (shape[k] <= nglue) and (procs[k] > 1):
                shape[k] *= 2
                procs[k] //= 2
                incr[k] *= 2
                r['g'+direc] = 1
                nglued += direc

        if len(nglued) == 1:
            # if there is one gluing, let's do another one
            nprocs = np.prod(procs)
            # print('only one glued, nprocs=%i' % nprocs)
            if nprocs == 1:
                # there is nothing we can do
                pass
            elif nprocs >= 2:
                direc = nglued
                k = 'kji'.index(direc)
                if procs[k] > 1:
                    shape[k] *= 2
                    procs[k] //= 2
                    incr[k] *= 2
                    r['g'+direc] = 2
                    # print('glue two instead of one')
                else:
                    # wait another turn
                    shape[k] //= 2
                    procs[k] *= 2
                    incr[k] //= 2
                    r.pop('g'+direc)
                    # print('postpone gluing to next level')

        grid = {'shape': shape.copy(), 'procs': procs.copy(),
                'incr': incr.copy(), 'coarsen': r.copy()}
        mygrids += [grid]
        lev += 1
        if np.prod(shape) <= ncellscoarsest:
            # assert np.prod(procs) == 1, 'weird problem with the grid hierachy'
            done = True

    return mygrids


def print_grids(grids):
    print('-'*80)
    print('Structure of the grids')
    print('======================')
    for lev, g in enumerate(grids):
        r = ''.join([s[1] for s in g['coarsen'].keys() if s[0] == 'r'])
        # glue = ''.join([s[1] for s in g['coarsen'].keys() if s[0] == 'g'])
        procs0 = [i*p for i, p in zip(g["incr"], g["procs"])]
        print('lev = %2i' % lev,
              '/ subdom size: %3i x %3i x %3i' % tuple(g['shape']),
              '/ cores: %1i x %1i x %1i' % tuple(g['procs']),
              '/ incr: %r' % (tuple(g["incr"]),),
              '/ procs0: %r' % (tuple(procs0),),
              '/ restrict: %3s' % r)
#              '/ glue: %3s' % glue)
    print('-'*80)


def set_mygluepartners(myrank, lev, grid: Grid, verbose=False):
    g = grid
    s = g.infos["subd"]
    # print("s=",s)
    glued = s["glued"]
    direc = s["gluing"]
    if grid.infos["glueflag"]:
        if verbose:
            print("** level %i" % lev, " / glue / procs= ", s["procs"])
        for d, r in glued.items():
            if type(r[0]) is list:
                r = [rr for rr in r if myrank in rr]
                if len(r) > 0:
                    r = r[0]
            if myrank in r:
                s = 2
                if (len(direc) == 1):
                    s = len(r)
                size = [1, 1, 1]
                for k, dd in enumerate('kji'):
                    if dd in direc:
                        size[k] = s

                mat = np.array(r).reshape(size)
                dom = d
                if verbose:
                    print("packed along %s :" %
                          direc, np.shape(mat), "\n", mat)
                k, j, i = np.where(mat == myrank)
                if verbose:
                    print("rank %i is located at " % myrank, k[0], j[0], i[0])
    else:
        if verbose:
            print("** level %i" % lev, " / no glue", s["procs"])
        mat = []
        dom = -1
    return mat, dom


class Dummygrid(Grid):
    def __init__(self, fine, coarse):
        assert any(coarse.incr != fine.incr)
        comm = MPI.COMM_WORLD
        nprocs = comm.Get_size()
        myrank = comm.Get_rank()
        matshape = tuple(coarse.incr // fine.incr)
        xshape = coarse.shape
        dummyshape = [j//i for i, j in zip(matshape, xshape)]
        ngbs = fine.neighbours
        nh = 1
        dummysize, domainindices = topo.get_variable_shape(
            dummyshape, ngbs, nh)
        N = np.prod(dummysize)

        self.shape = dummyshape
        self.neighbours = ngbs
        self.size = dummysize
        self.N = N
        self.domainindices = domainindices
        self.nh = nh
        self.x = np.zeros((self.N,))
        self.A = []

    def toarray(self, which):
        """
        Return the 3D array view of 'x'
        """
        assert which == "x"
        return np.reshape(getattr(self, which), self.size)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    import topology as topo

    procs = [8, 4, 2]
    topology = 'closed'
    topo.topology = topology

    param = {'nx': 8, 'ny': 128, 'nz': 64, 'nh': 1,
             'procs': procs, 'nglue': 16, 'ncellscoarsest': 16}

    allgrids = define_grids_v2(param)
    print_grids(allgrids)
