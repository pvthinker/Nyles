"""

set of functions related to subdomain partition

  - where I am in the global domain?
  - who are my neighbours?
  - what is the width of my halo?
  - ...

"""

import numpy as np
import itertools as iter

# this is a global variable
# assuming this module is imported as
# > import topology as topo
# it is accessed via
# > topo.topology
# to change its value in all modules that
# import it, one simply redefine
# topo.topology = 'closed'
# at only one place and all modules
# will see the change

topology = 'closed'


def rank2loc(rank, procs):
    """ Return the tuple (k, j, i) of subdomain location

    - rank, the MPI rank
    - procs, the tuple of number of subdomains in each direction

    """

    loc = [rank // (procs[2]*procs[1]),
           (rank // procs[2]) % procs[1],
           rank % procs[2]]

    return loc


def loc2rank(loc, procs):
    """ Return the MPI rank using

    - loc, the tuple of subdomain location
    - procs, the tuple of number of subdomains in each direction

    loc2rank is the reverse operation of rank2loc

    """
    rank = (loc[0]*procs[1] + loc[1])*procs[2] + loc[2]
    return rank

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
    return loc2rank([k, j, i], proc)


def get_neighbours(location, procs, incr=[1, 1, 1], extension=26):
    """ Return the neighbours rank in the six directions

    the result is presented in a dictionnary

    topology = 'closed' or 'perio_x', 'perio_xy', 'perio_y', 'perio_xyz' ...

    None stands for no neighbour


    """
    if type(location) is int:
        k, j, i = rank2loc(location, procs)
    elif type(location) is list:
        k, j, i = location
    else:
        raise ValueError
    
    nz, ny, nx = procs
    incz, incy, incx = incr
    k *= incz
    j *= incy
    i *= incx

    # alldirec has 27 elements, the coordinates of
    # the 3x3 cube
    alldirec = [(a, b, c) for a, b, c in iter.product(
        [-1, 0, 1], [-1, 0, 1], [-1, 0, 1])]

    # depending on 'extension' we extract the
    # appropriate directions
    if extension == 6:
        directions = [d for d in alldirec if sum(np.abs(d)) == 1]

    elif extension == 18:
        directions = [d for d in alldirec if sum(np.abs(d)) in [1, 2]]

    elif extension == 26:
        directions = [d for d in alldirec if sum(np.abs(d)) >= 1]

    else:
        raise ValueError('neighbours should be 6, 18 or 26')

    ngs = {}
    for dk, dj, di in directions:
        ng = loc2rank([(k+dk*incz) % nz,
                       (j+dj*incy) % ny,
                       (i+di*incx) % nx], procs)
        if 'x' in topology:
            pass
        else:
            if (i+di*incx) < 0 or (i+di*incx) >= nx:
                ng = None

        if 'y' in topology:
            pass
        else:
            if (j+dj*incy) < 0 or (j+dj*incy) >= ny:
                ng = None

        if 'z' in topology:
            pass
        else:
            if (k+dk*incz) < 0 or (k+dk*incy) >= nz:
                ng = None

        if ng is None:
            # don't even keep track of that direction
            pass
        else:
            # neighbour key is the tuple of direction
            # not a plain string because with up to
            # 28 neighbours naming all of them is rather
            # cumbersome
            ngs[(dk, dj, di)] = ng
    return ngs


def get_variable_shape(innersize, ngbs, nh):
    """ 
    innersize sets the interior domain size
    innersize = [nz, ny, nx]

    ngbs is the neighbours dictionary

    return:

    size = the extended domain size (with halo)

    domainindices = (k0, k1, j0, j1, i0, i1),
                    the list of start and last index
    """
    size = innersize.copy()

    if (-1, 0, 0) in ngbs.keys():
        size[0] += nh
        k0 = nh
    else:
        size[0] += 1
        k0 = 1
    k1 = size[0]
    if (+1, 0, 0) in ngbs.keys():
        size[0] += nh

    if (0, -1, 0) in ngbs.keys():
        size[1] += nh
        j0 = nh
    else:
        size[1] += 1
        j0 = 1
    j1 = size[1]
    if (0, +1, 0) in ngbs.keys():
        size[1] += nh

    if (0, 0, -1) in ngbs.keys():
        size[2] += nh
        i0 = nh
    else:
        size[2] += 1
        i0 = 1
    i1 = size[2]
    if (0, 0, +1) in ngbs.keys():
        size[2] += nh

    # k0 is the vertical index of the first interior point
    # in python indexing convention (zero starting).
    # points with k<k0 are in the halo
    # if k0 == 0, there is no halo on this side
    # k1 is the vertical index of the last inzlnterior point
    # points with k>k1 are in the halo
    # if k1 = nzl-1, there is no halo on this side
    domainindices = (k0, k1, j0, j1, i0, i1)
    return size, domainindices

def check_graph(allneighbours):
    """
    build the connectivity matrix of the neighbours graph

    assert that it is symmetric

    """
    nbprocs = len(allneighbours)
    assert all([type(n) == dict for n in allneighbours]
               ), "connectivity can be checked only with the list of neighbours for *all* ranks"

    connectivity = np.zeros((nbprocs, nbprocs), dtype=int)
    for i, n in enumerate(allneighbours):
        print(n)
        for key, j in n.items():
            if j is None:
                pass
            else:
                connectivity[i, j] += 1

    print(connectivity)

    # the matrix should be symmetric
    msg = 'neighbours are not connected symmetrically'
    assert (connectivity == connectivity.transpose()).all(), msg
    print('connectivity matrix is symmetric')

def get_variable_shape(innersize, ngbs, nh):
    """ 
    innersize sets the interior domain size
    innersize = [nz, ny, nx]

    ngbs is the neighbours dictionary

    return:

    size = the extended domain size (with halo)

    domainindices = (k0, k1, j0, j1, i0, i1),
                    the list of start and last index
    """
    size = innersize.copy()

    if (-1, 0, 0) in ngbs.keys():
        size[0] += nh
        k0 = nh
    else:
        k0 = 0
    k1 = size[0]
    if (+1, 0, 0) in ngbs.keys():
        size[0] += nh

    if (0, -1, 0) in ngbs.keys():
        size[1] += nh
        j0 = nh
    else:
        j0 = 0
    j1 = size[1]
    if (0, +1, 0) in ngbs.keys():
        size[1] += nh

    if (0, 0, -1) in ngbs.keys():
        size[2] += nh
        i0 = nh
    else:
        i0 = 0
    i1 = size[2]
    if (0, 0, +1) in ngbs.keys():
        size[2] += nh

    # k0 is the vertical index of the first interior point
    # in python indexing convention (zero starting).
    # points with k<k0 are in the halo
    # if k0 == 0, there is no halo on this side
    # k1 is the vertical index of the last inzlnterior point
    # points with k>k1 are in the halo
    # if k1 = nzl-1, there is no halo on this side
    domainindices = (k0, k1, j0, j1, i0, i1)
    return size, domainindices


if __name__ == '__main__':

    # procs = [4, 2, 1] stands for 4 subdomains in z, 2 in y, 1 in x

    procs = [1, 4, 4]
    topology = 'perio_x'

    ngs = []
    for myrank in range(np.prod(procs)):
        loc = rank2loc(myrank, procs)
        neighbours = get_neighbours(loc, procs, extension=18)

        print('-'*60)
        print('rank %i, loc (%i, %i, %i)' % (myrank, *loc))
        print(' neighbours: ')
        for k, v in neighbours.items():
            print('    ', k, v)
        ngs += [neighbours]
    print('-'*60)
    check_graph(ngs)
