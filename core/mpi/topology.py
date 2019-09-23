"""

set of functions related to subdomain partition

  - where I am in the global domain?
  - who are my neighbours?
  - what is the width of my halo?
  - ...

"""

import numpy as np


def rank2loc(rank, procs):
    """ Return the tuple (k, j, i) of subdomain location

    - rank, the MPI rank
    - procs, the tuple of number of subdomains in each direction

    """

    nbproc = np.prod(procs)
    ranks = np.reshape(np.arange(nbproc), procs)
    loc = [idx[0] for idx in np.where(ranks == rank)]
    return loc


def loc2rank(loc, procs):
    """ Return the MPI rank using

    - loc, the tuple of subdomain location
    - procs, the tuple of number of subdomains in each direction

    loc2rank is the reverse operation of rank2loc

    """
    nbproc = np.prod(procs)
    ranks = np.reshape(np.arange(nbproc), procs)
    rank = ranks[loc[0], loc[1], loc[2]]
    return rank


def get_neighbours(loc, procs, topology):
    """ Return the neighbours rank in the six directions

    the result is presented in a dictionnary

    topology = 'closed' or 'perio_x', 'perio_xy', 'perio_y', 'perio_xyz' ...

    None stands for no neighbour


    """
    k, j, i = loc
    nz, ny, nx = procs

    im = loc2rank([k, j, (i-1) % nx], procs)
    ip = loc2rank([k, j, (i+1) % nx], procs)
    jm = loc2rank([k, (j-1) % ny, i], procs)
    jp = loc2rank([k, (j+1) % ny, i], procs)
    km = loc2rank([(k-1) % nz, j, i], procs)
    kp = loc2rank([(k+1) % nz, j, i], procs)

    if 'x' in topology:
        pass
    else:
        if i == 0:
            im = None
        if i == nx-1:
            ip = None

    if 'y' in topology:
        pass
    else:
        if j == 0:
            jm = None
        if j == ny-1:
            jp = None

    if 'z' in topology:
        pass
    else:
        if k == 0:
            km = None
        if k == nz-1:
            kp = None

    return {'km': km, 'kp': kp, 'jm': jm, 'jp': jp, 'im': im, 'ip': ip}


def get_halowidth(myrank, procs, topology, nh=2):
    """
    return the halo width on the left and the right for each direction
    in the form of a dictionary

    width = 0, if no neightbour ! THIS IS VERY NON CONVENTIONAL
    width = nh, otherwise.

    """
    loc = rank2loc(myrank, procs)
    neighbours = get_neighbours(loc, procs, topology)
    halowidth = {}
    for key, rank in neighbours.items():
        if rank is None:
            halowidth[key] = 0
        else:
            halowidth[key] = nh
    return halowidth


def noneighbours():
    """
    return a dictionnary with no neighbour in all directions
    used in the monoprocessor case, topology = 'closed'
    """
    neighbours = {}
    for k in ['km', 'kp', 'jm', 'jp', 'im', 'ip']:
        neighbours[k] = None
    return neighbours


if __name__ == '__main__':
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nbproc = MPI.COMM_WORLD.Get_size()

    nh = 3
    nx = 64
    ny = 48
    nz = 30

    # procs = [4, 2, 1] stands for 4 subdomains in z, 2 in y, 1 in x
    # consequently the code should be called with 8 processes
    procs = [4, 2, 1]
    topology = 'perio_y'

    loc = rank2loc(myrank, procs)
    neighbours = get_neighbours(loc, procs, topology)

    halowidth = get_halowidth(myrank, procs, topology, nh=nh)

    print('-'*60)
    print('rank %i, loc (%i, %i, %i)' % (myrank, *loc))
    print(' halowidth:  ', halowidth)
    print(' neighbours: ', neighbours)
