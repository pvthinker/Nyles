"""

Tools that implement the halo filling

"""

import itertools as itert
import numpy as np
from mpi4py import MPI
import topology as topo
#import fortran_halo as fortran
from timing import timing


class Halo():
    def __init__(self, grid):
        """

        scalar is an instance of the Scalar object

        create a halo method once we have a scalar available
        e.g. 'b' the buoyancy

        b = state.get('b')
        halo = Halo(b)

        we can now fill any Scalar or Vector with

        halo.fill(b)
        halo.fill(u)

        where u is a vector, e.g. u = state.get('u')

        halo defines it own preallocated buffers and MPI requests

        The local topology is known from param['neighbours'] stored
        in the scalar. It is worth saying that each halo only needs
        the local information = who are my neighbours. A halo does
        not need to know the global topology.

        """
        comm = MPI.COMM_WORLD
        myrank = comm.Get_rank()
        neighbours = grid['neighbours']
        nh = grid['nh']
        #extension = grid['extension']
        self.size = grid['size']
        self.domainindices = grid['domainindices']
        self.neighbours = neighbours
        self.myrank = myrank
        self.nh = nh
        nz, ny, nx = grid['shape']
        shape = grid['shape']
        self.nz, self.ny, self.nx = nz, ny, nx

        #icase = {6: 1, 18: 2, 26: 3}
        #assert extension in icase.keys()
        #self.icase = icase[extension]

        alldirec = [(a, b, c) for a, b, c in itert.product(
            [-1, 0, 1], [-1, 0, 1], [-1, 0, 1]) if abs(a)+abs(b)+abs(c) > 0]

        # allocate buffers
        self.sbuf = {}
        self.rbuf = {}
        for direc in neighbours.keys():
            bufshape = [nh, nh, nh]
            for l in range(3):
                if abs(direc[l]) == 0:
                    bufshape[l] = shape[l]
            self.sbuf[direc] = np.zeros(bufshape)
            self.rbuf[direc] = np.zeros(bufshape)

        # define MPI requests
        self.reqs = []
        self.reqr = []
        tag = 0
        for direc, yourrank in neighbours.items():
            # reqs = comm.Send_init(self.sbuf[direc], yourrank, tag=myrank)
            # reqr = comm.Recv_init(self.rbuf[direc], yourrank, tag=yourrank)
            flipdirec = tuple([-k for k in direc])
            tag = 0
            ftag = 0
            for k, d in enumerate(direc):
                tag += (3**k)*(d+1)
                ftag += (3**k)*(-d+1)
            reqs = comm.Send_init(self.sbuf[direc], yourrank, tag=myrank+tag)
            reqr = comm.Recv_init(self.rbuf[direc], yourrank, tag=yourrank+ftag)
            self.reqs += [reqs]
            self.reqr += [reqr]
            tag += 1

        # define interior domain slices
        domi = self.domainindices
        self.iidx = []
        # self.iidx is a 3 x 3 combined list x dictionary of slices
        # the list index refers to the dimension 0: k, 1: j, 2: i
        # the dictionnary entry refers to the direction -1, 0 or 1
        # slice sweeps in the interior domain = places where to grab
        # known x to be transfered in the neighbour's halo
        for l in range(3):
            idx = {}
            idx[-1] = slice(-nh-nh, -nh)
            idx[1] = slice(nh, nh+nh)
            idx[0] = slice(domi[2*l], domi[2*l+1])
            self.iidx += [idx]

        # define outer domain slices
        domi = self.domainindices
        size = self.size
        self.oidx = []
        # self.oidx has the same structure
        # slice sweeps in the halo domain = places where to fill halo
        # with known x from the neighbour's interior domain
        for l in range(3):
            idx = {}
            idx[-1] = slice(0, nh)
            idx[1] = slice(self.size[l]-nh, self.size[l])
            idx[0] = slice(domi[2*l], domi[2*l+1])
            self.oidx += [idx]

#    @timing
    def fill(self, thing):
        """
        thing is either a scalar, a numpy array or a vector
        """
        nature = type(thing).__name__
        if nature == 'Scalar':
            self.fillarray(thing.view('i'))

        elif nature == 'ndarray':
            self.fillarray(thing)

        elif nature == 'Vector':
            self.fillvector(thing)

        else:
            raise ValueError('try to fill halo with unidentified object')

    def fillarray(self, x):
        """
        this is where all the MPI instructions are

        there are basically four steps:

        - 1) copy inner values of x into buffers
        - 2) send buffers
        - 3) receive buffers
        - 4) copy buffers into x halo

        we use predefined requests and lauch them with
        Prequest.Startall(). This must be used in
        conjunction with Prequest.Waitall()
        that ensures all requests have been completed

        """

        MPI.Prequest.Startall(self.reqr)

        # 1) halo to buffer
        idx = self.iidx
        for direc in self.neighbours.keys():
            dk, dj, di = direc
            self.sbuf[direc][:, :, :] = x[idx[0][-dk], idx[1][-dj], idx[2][-di]]

        # 2)
        MPI.Prequest.Startall(self.reqs)

        # 3)
        MPI.Prequest.Waitall(self.reqr)

        # 4) buffer to halo
        idx = self.oidx
        for direc in self.neighbours.keys():
            dk, dj, di = direc
            x[idx[0][dk], idx[1][dj], idx[2][di]] = self.rbuf[direc]

        MPI.Prequest.Waitall(self.reqs)

    def fillvector(self, vector):
        """
        trivial extension to a vector

        TODO: see if it's worth optimizing
        """
        for direc in 'ijk':
            u = vector[direc].view('i')
            self.fillarray(u)

    def _print(self):
        """
        Print the shape of send buffers

        (for debug purpose)
        """
        for key, buf in self.sbuf.items():
            print(self.myrank, key, np.shape(buf))

def set_halo(param, state):
    extension = 6
    nh = param["nh"]
    procs = param["procs"]
    neighbours = param["neighbours"]
    shape = state.b.shape
    size, domainindices = np.shape(state.b.view('i')), state.b.domainindices
    localgrid = {'shape': shape,
        'size': size,
        'nh': nh,
        'neighbours': neighbours,
        'domainindices': domainindices,
        'extension': extension}
    check_halo_width(procs, shape, nh)
    return Halo(localgrid)

def test_111_domain():
    topo.topology = 'perio_xyz'
    extension = 26
    nz = 1
    ny = 1
    nx = 1
    nh = 1

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs, extension=extension)

    param = {'nx': nx, 'ny': ny, 'nz': nz, 'nh': nh, 'neighbours': neighbours}

    state = var.get_state(param)

    b = state.get('b')
    # we can define the halo toolbox

    grid = {'shape': [nz, ny, nx],
            'size': b.size,
            'nh': nh,
            'neighbours': neighbours,
            'domainindices': b.domainindices,
            'extension': extension}

    halo = Halo(grid)

    # for r, buf in halo.rbuf.items():
    #     if r in neighbours.keys():
    #         buf[:] = neighbours[r]*1.

    #     else:
    #         buf[:] = -1.

    # set the field equal to myrank
    # making clear that the halo is wrong
    b.activeview = 'i'
    # if myrank == 0:
    #     b.view('i')[:, :, :] = 999
    # else:
    #     b.view('i')[:, :, :] = 0.
    x = b.view('i')

    x[:, :, :] = myrank

    halo.fill(b)
    # now in the halo we'd have the rank of the neighbour

    y = np.asarray(x.copy(), dtype=int)
    msg = 'pb withneighbours'
    for direc, yourrank in neighbours.items():
        dk, dj, di = direc
        assert(y[1+dk, 1+dj, 1+di] == yourrank), msg
    if myrank == 0:
        print('fill halo of Scalar is okay')

    l = nz+2*nh
    m = ny+2*nh
    n = nx+2*nh
    x = np.zeros((l, m, n))

    grid = {'shape': [nz, ny, nx],
            'size': np.shape(x),
            'nh': nh,
            'neighbours': neighbours,
            'domainindices': b.domainindices,
            'extension': extension}

    halox = Halo(grid)
    x[:] = myrank
    halox.fill(x)
    y = np.asarray(x.copy(), dtype=int)
    for direc, yourrank in neighbours.items():
        dk, dj, di = direc
        #print(myrank, y[1+dk, 1+dj, 1+di] , yourrank)
        assert(y[1+dk, 1+dj, 1+di] == yourrank), msg
    if myrank == 0:
        print('fill halo of ndarray is okay')


def check_halo_width(procs, shape, nh):
    msg = 'domain is too narrow to cope with a halo that wide'
    if (procs[0] > 1) or ('z' in topo.topology):
        assert shape[0] >= nh, msg
    if (procs[1] > 1) or ('y' in topo.topology):
        assert shape[1] >= nh, msg
    if (procs[2] > 1) or ('x' in topo.topology):
        msg2 = msg + ' / shape[2]={} < nh'.format(shape[2])
        assert shape[2] >= nh, msg2


if __name__ == '__main__':

    import variables as var
    import sys

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    procs = [1, 3, 3]

    msg = 'use mpirun -np %i python ' % np.prod(procs) + ' '.join(sys.argv)
    assert comm.Get_size() == np.prod(procs), msg

    #test_111_domain()


    topo.topology = 'closed'
    extension = 18
    nz, ny, nx = 1, 2, 2
    shape = [nz, ny, nx]
    nh = 2
    loc = topo.rank2loc(myrank, procs)
    print('myrank=%i / loc=' % myrank, loc)
    neighbours = topo.get_neighbours(loc, procs, extension=extension)
    size, domainindices = topo.get_variable_shape([nz, ny, nx], neighbours, nh)
    grid = {'shape': [nz, ny, nx],
            'size': size,
            'nh': nh,
            'neighbours': neighbours,
            'domainindices': domainindices,
            'extension': extension}
    check_halo_width(procs, shape, nh)
    halox = Halo(grid)
    x = np.zeros(size)
    x[:] = myrank
    halox.fill(x)

    if myrank in [4, 2]:
        for k in range(np.shape(x)[0]):
            print('rank=%i / k=%i' % (myrank, k), domainindices)
            print(x[k])
