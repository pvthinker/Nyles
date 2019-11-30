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

class Grid(object):
    """

    Provide functions working on a particular grid

      - norm
      - smooth
      - residual
      - solveexactly

    """
    def __init__(self, grid, param):
        for key in ['shape', 'neighbours', 'nh', 'size', 'N', 'domainindices']:
            setattr(self, key, grid[key])
            
        for key in ['omega', 'ndeepest']:
            setattr(self, key, param[key])

        self.x = np.zeros((self.N,))
        self.b = np.zeros((self.N,))
        self.r = np.zeros((self.N,))
        self.halo = halo.Halo(grid)

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

    def norm(self, which='r'):
        """
        Compute the norm of 'which'
        excluding the halo
        """
        y = self.toarray(which)
        k0, k1, j0, j1, i0, i1 = self.domainindices
        y2 = fortran.norm(y, k0, k1, j0, j1, i0, i1)
        # todo : add a global communication to sum y2
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
    nh = param['nh']

    # todo: 'nglue' and 'ncellscoarsest' should be set from 'param'
    nglue = 16
    ncellscoarsest = 32

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
            assert np.prod(procs) == 1, 'weird problem with the grid hierachy'
            done = True

    return mygrids


def print_grids(grids):
    print('-'*80)
    print('Structure of the grids')
    print('======================')
    for lev, g in enumerate(grids):
        r = ''.join([s[1] for s in g['coarsen'].keys() if s[0] == 'r'])
        glue = ''.join([s[1] for s in g['coarsen'].keys() if s[0] == 'g'])
        print('lev = %2i' % lev,
              '/ subdom size: %3i x %3i x %3i' % tuple(g['shape']),
              '/ cores: %1i x %1i x %1i' % tuple(g['procs']),
              '/ restrict: %3s' % r,
              '/ glue: %3s' % glue)
    print('-'*80)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    import topology as topo

    procs = [2, 2, 2]
    topology = 'closed'
    myrank = 3
    nh = 3
    topology = 'closed'

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    param = {'nx': 128, 'ny': 128, 'nz': 64, 'nh': nh,
             'neighbours': neighbours, 'procs': procs}
