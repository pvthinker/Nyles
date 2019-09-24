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

"""
import numpy as np
import topology as topo


def define_grids(p):
    """
    define the grid hierarchy along with the coarsening/gluing operations to do

    gluing is triggered when the size in that direction is <= 'nglue'

    a grid is declared the coarsest when all cores compute the whole domain
    and that the number of cells of the global domain is <= 'ncellscoarsest'
    """
    procs = p['procs']
    nx, ny, nz = p['nx'], p['ny'], p['nz']
    shape = [nz, ny, nx]

    nglue = 16
    ncellscoarsest = 16

    lev = 0
    coarse = []
    done = False
    grids = [shape.copy()]
    cores = [procs.copy()]
    while not(done):
        r = {}
        done = True
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

        grids += [shape.copy()]
        coarse += [r]
        cores += [procs.copy()]
        lev += 1
        if np.prod(shape) <= ncellscoarsest:
            assert np.prod(procs) == 1, 'weird problem with the grid hierachy'
            done = True

    nblevs = lev+1
    print('nb levels %i' % nblevs)
    return grids, coarse, cores


def ranktoloc(rank, core):
    pk, pj, pi = core
    i = (rank) % pi
    j = (rank//pi) % pj
    k = (rank//(pi*pj)) % pk
    return [k, j, i]


def myloc(loc, inc, delta, proc):
    """
    generalized loc2rank with incr != 1

    used to determine the neighbours

    """
    k, j, i = loc
    incz, incy, incx = inc
    k = (k+delta[0])*incz % proc[0]
    j = (j+delta[1])*incy % proc[1]
    i = (i+delta[2])*incx % proc[2]
    return topo.loc2rank([k, j, i], proc)


def pairing(coarse, cores):
    """

    describe the gluing process and the subdomains on the intermediate grids

    """
    procs = cores[0].copy()
    procs0 = procs.copy()
    nb = np.prod(procs)
    incr = [1, 1, 1]
    rank = np.arange(nb)
    loc = ranktoloc(rank, procs)
    loc0 = loc.copy()
    k0, j0, i0 = loc0
    gathers = []
    doms = []
    flag = True
    first = True
    tags = rank*0

    print('-'*80)
    print('')
    print('Description of each subdomain decomposition')
    print('===========================================')
    print('')
    print('starting from the finest level: one rank per subdomain')
    print('  ending at the coarsest level: all rank on the whole domain')
    print('')
    for lev, c in enumerate(coarse):
        keys = c.keys()
        if flag:
            loc0 = loc.copy()
            k0, j0, i0 = loc0
        if any(['g' in k for k in keys]) or lev == 0:
            print('-'*80)
            g = []
            for n, direc in enumerate('kji'):
                if 'g'+direc in keys:
                    incr[n] *= 2
                    procs[n] //= 2
                    loc[n] //= 2
                    g += [direc]
            if g != '':
                gathers += [g]
                flag = True
            else:
                flag = False
                print('***no gather****')
            k, j, i = loc

            nz, ny, nx = procs0
            incz, incy, incx = incr
            dom = [0]*nb
            for r in rank:
                l = [k[r], j[r], i[r]]
                # d is the index of the subdomain
                # rank with same d do the same computation
                d = myloc(l, incr, [0, 0, 0], procs0)

                # here is the list of neighbours
                #
                # the print has been discarded to lighten the output
                # but clearly this information is very important
                # note that the neighbours relationships for a given
                # rank depends on the subdomain decomposition
                # the neighbours get further and further as
                # subdomains are glued. This is accounted by
                # 'incr'.  incr == 1 on the finest grid, then
                # whenever there is gluing in one direction
                # incr is doubled in that direction, until
                # the point where the is only one subdomain is that
                # direction. Each rank is then neighbour with itself
                #
                im = myloc(l, incr, [0, 0, -1], procs0)
                ip = myloc(l, incr, [0, 0, +1], procs0)
                jm = myloc(l, incr, [0, -1, 0], procs0)
                jp = myloc(l, incr, [0, +1, 0], procs0)
                km = myloc(l, incr, [-1, 0, 0], procs0)
                kp = myloc(l, incr, [+1, 0, 0], procs0)
                dom[r] = d

                # discarded print
                # print('r: %i - d: %i / ' % (r, d), im, ip, jm, jp, km, kp)

            doms += [dom]
            dd = set(dom)
            print('Subdomains are : ', dd)
            nbproc = np.prod(procs)
            print('corresponding to domain decomposition: ', procs)
            print('increments in each directions: ', incr)

            print('Location of each subdomain')
            print('k: ', k)
            print('j: ', j)
            print('i: ', i)
            if not(first):
                print('gluing was in directions ', g)
            for d in dd:
                family = [r for r in rank if dom[r] == d]
                # in a family, all ranks have the same loc
                # because they all compute the same subdomain
                locs = [(k0[r], j0[r], i0[r]) for r in family]

                # tags is defined from the previous domain decomposition
                glued = [[r for r in family if tags[r] == t]
                         for t in set(tags[family])]

                print(' - subdomain : ', set([d]))
                print('     located at ', locs[0])
                print('     computed by ' % d, family)
                if not(first):
                    print('     obtained by gluing ranks:', glued)
            #
            # tags is used to identify which ranks are glued together
            #
            # glued ranks had identical tag in the previous
            # subdomain decomposition
            #
            for d in dd:
                family = [r for r in rank if dom[r] == d]
                for k, r in enumerate(family):
                    tags[r] = k
            first = False
    print('-'*80)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    procs = [2, 4, 2]
    topology = 'closed'
    myrank = 3
    nh = 3

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs, topology)

    param = {'nx': 64, 'ny': 128, 'nz': 64, 'nh': nh,
             'neighbours': neighbours, 'procs': procs}

    print('-'*80)
    grids, coarse, cores = define_grids(param)
    print('-'*80)
    print('Structure of the grids')
    for lev, (g, c) in enumerate(zip(grids, cores)):
        print('lev = %2i' % lev, 'subdom size: ', g, 'subdom decomp: ', c)
    print('-'*80)
    print('Coarsening operations')
    for k, c in enumerate(coarse):
        print('from lev = %i to %i' % (k, k+1),
              ' coarsening in ', [k[1] for k in c.keys() if k[0] == 'r'],
              ' gluing in ', [k[1] for k in c.keys() if k[0] == 'g'])

    pairing(coarse, cores)
