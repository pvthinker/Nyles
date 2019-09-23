"""

Tools that implement the halo filling

"""

from mpi4py import MPI
import numpy as np
from mpi import topology as topo


class Halo(object):
    def __init__(self, scalar):
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
        neighbours = scalar.param['neighbours']
        nh = scalar.param['nh']
        size = scalar.size

        self.myrank = myrank
        self.nh = nh

        self.sbuf = {}
        self.rbuf = {}
        self.reqs = {}
        self.reqr = {}
        self.slab = {}
        for direc in neighbours.keys():
            yourrank = neighbours[direc]
            if yourrank is None:
                slab = [0, 0]

            else:
                if direc[0] == 'i':
                    slab = [size[2], size[1]]

                if direc[0] == 'j':
                    slab = [size[0], size[2]]

                if direc[0] == 'k':
                    slab = [size[1], size[0]]

                shape = [nh] + slab
                sbuf = np.zeros(shape)
                rbuf = np.zeros(shape)
                reqs = comm.Send_init(sbuf, yourrank, myrank)
                reqr = comm.Recv_init(rbuf, yourrank, yourrank)

                self.sbuf[direc] = sbuf
                self.rbuf[direc] = rbuf
                self.reqs[direc] = reqs
                self.reqr[direc] = reqr

            self.slab[direc] = slab

    def fill(self, thing):
        """

        thing is either a scalar or a vector

        """
        nature = type(thing).__name__
        if nature == 'Scalar':
            self.fillscalar(thing)

        elif nature == 'Vector':
            self.fillvector(thing)

        else:
            raise ValueError('try to fill halo with unidentified object')

    def fillscalar(self, scalar):
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
        nh = self.nh

        # convert dictionary into list
        reqs = [req for key, req in self.reqs.items()]
        reqr = [req for key, req in self.reqr.items()]

        MPI.Prequest.Startall(reqr)
        # 1) halo to buffer
        for direc in 'ijk':
            x = scalar.flipview(direc)
            sm = np.prod(self.slab[direc+'m'])
            sp = np.prod(self.slab[direc+'p'])
            # check whether there is a left and/or a right halo
            if (sm > 0) & (sp > 0):
                bm = self.sbuf[direc+'m']
                bp = self.sbuf[direc+'p']
                bm[:, :, :] = x[nh:2*nh, :, :]
                bp[:, :, :] = x[-2*nh:-nh, :, :]

            elif (sm > 0) & (sp == 0):
                bm = self.sbuf[direc+'m']
                bm[:, :, :] = x[nh:2*nh, :, :]

            elif (sm == 0) & (sp > 0):
                bp = self.sbuf[direc+'p']
                bp[:, :, :] = x[-2*nh:-nh, :, :]

            else:
                pass
        # 2)
        MPI.Prequest.Startall(reqs)

        # 3)
        MPI.Prequest.Waitall(reqr)
        # ierr = 9999
        # try:
        #     ierr = MPI.Prequest.Waitall(reqr)
        # except:
        #     print('IERR', ierr)

        # 4) buffer to halo
        for direc in 'ijk':
            x = scalar.flipview(direc)
            sm = np.prod(self.slab[direc+'m'])
            sp = np.prod(self.slab[direc+'p'])
            if (sm > 0) & (sp > 0):
                bm = self.rbuf[direc+'m']
                bp = self.rbuf[direc+'p']
                x[:nh, :, :] = bm
                x[-nh:, :, :] = bp

            elif (sm > 0) & (sp == 0):
                bm = self.rbuf[direc+'m']
                x[:nh, :, :] = bm

            elif (sm == 0) & (sp > 0):
                bp = self.rbuf[direc+'p']
                x[-nh:, :, :] = bp

            else:
                pass

        MPI.Prequest.Waitall(reqs)

    def fillvector(self, vector):
        """
        trivial extension to a vector

        TODO: see if it's worth optimizing
        """
        for direc in 'ijk':
            self.fillscalar(getattr(vector, direc))

    def _print(self):
        """
        Print the shape of send buffers

        (for debug purpose)
        """
        for key, buf in self.sbuf.items():
            print(self.myrank, key, np.shape(buf))


if __name__ == '__main__':

    import variables as var

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    procs = [1, 2, 2]
    topology = 'perio_y'
    nh = 3

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs, topology)

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': nh,
             'neighbours': neighbours}

    state = var.get_state(param)

    b = state.get('b')
    # we can define the halo toolbox
    halo = Halo(b)

    # set the field equal to myrank
    # making clear that the halo is wrong
    b.view('i')[:, :, :] = myrank
    halo.fill(b)
    # now in the halo we'd have the rank of the neighbour
    x = b.view('i')
    print('myrank %i' % myrank, x[5, 5, :])

    # check that halo.fill() also accepts a vector
    u = state.get('u')
    halo.fill(u)
