""" 

Tools that implement the halo filling

"""

from mpi4py import MPI
import numpy as np
import timeit


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


def get_neighbours(loc, procs, topology='cube'):
    """ Return the neighbours rank in the six directions

    the result is presented in a dictionnary

    None stands for no neighbour

    """
    k, j, i = loc
    nz, ny, nx = procs

    if topology == 'cube':
        im = loc2rank([k, j, (i-1) % nx], procs)
        ip = loc2rank([k, j, (i+1) % nx], procs)
        jm = loc2rank([k, (j-1) % ny, i], procs)
        jp = loc2rank([k, (j+1) % ny, i], procs)
        km = loc2rank([(k-1) % nz, j, i], procs)
        kp = loc2rank([(k+1) % nz, j, i], procs)
        if nx == 1:
            im, ip = None, None

        if ny == 1:
            jm, jp = None, None

        if k == 0:
            km = None
        if k == nz-1:
            kp = None
    return {'km': km, 'kp': kp, 'jm': jm, 'jp': jp, 'im': im, 'ip': ip}

class Halo_outerdim(object):
    # STATUS: it's a try, the original one is Halo(object)
    # Idea behind: do the halo only on the two hyperplanes
    # (the first and the last, within the convention)
    """ 
    Class to fill the halo only in the outer dimension

    The halo to update is therefore contiguous in memory

    the bottom halo is x[:nh, :, :]
    the top halo is x[-nh:, :, :]

    they get their values from respectively
    x[nh:nh+nh, :, :] and x[-nh-nh:-nh, :, :]

    - nh is the halo width

    """
    def __init__(self, comm, grid, neighbours, outerdir='k'):
        self.nh = grid['nh']
        self.shape = grid['shape']
        self.procs = grid['procs']
        self.comm = comm
        self.neighbours = neighbours

        myrank = comm.Get_rank()
        self.myrank = myrank

        nz, ny, nx = self.shape
        if outerdir == 'k':
            n2, n1, n0 = nz, ny, nx
            
        elif outerdir == 'j':
            n2, n1, n0 = ny, nx, nz
            
        elif outerdir == 'i':
            n2, n1, n0 = nx, nz, ny
            
        nh = self.nh
        nghb = neighbours
        
        self.sbuf = []
        self.rbuf = []
        self.reqs = []
        self.reqr = []
        self.nghbs = []
        self.direcs = []
        for direc in neighbours.keys():
            yourrank = neighbours[direc]
            if (yourrank is None) or direc[0] != outerdir:
                pass
            else:
                shape = (nh, n1, n0)
                sbuf = np.zeros(shape)
                rbuf = np.zeros(shape)
                reqs = comm.Send_init(sbuf, yourrank, myrank)
                reqr = comm.Recv_init(rbuf, yourrank, yourrank)

                self.sbuf += [sbuf]
                self.rbuf += [rbuf]
                self.reqs += [reqs]
                self.reqr += [reqr]
                self.nghbs += [yourrank]
                self.direcs += [direc]

    def fill(self, x):
        """ fill the halo of x"""

        reqs = self.reqs
        reqr = self.reqr

        MPI.Prequest.Startall(reqr)
        # halo to buffer
        for buf, direc in zip(self.sbuf, self.direcs):
            if direc[1] == 'm':
                buf[:, :] = x[nh:nh+nh, :, :]
            if direc[1] == 'p':
                buf[:, :] = x[-nh-nh:-nh, :, :]
            
        MPI.Prequest.Startall(reqs)
        ierr = 9999
        try:
            ierr = MPI.Prequest.Waitall(reqr)
        except:
            print('IERR', ierr)
        # buffer to halo
        for buf, direc in zip(self.rbuf, self.direcs):
            if direc[1] == 'm':
                x[:nh, :, :] = buf[:, :]
            if direc[1] == 'p':
                x[-nh:, :, :] = buf[:, :]

        MPI.Prequest.Waitall(reqs)
        
    
class Halo(object):
    def __init__(self, comm, grid, neighbours):
        self.nh = grid['nh']
        self.shape = grid['shape']
        self.procs = grid['procs']
        self.comm = comm
        self.neighbours = neighbours

        myrank = comm.Get_rank()
        self.myrank = myrank

        nz, ny, nx = self.shape
        nh = self.nh
        nghb = neighbours

        zface = (nx-2*nh)*(ny-2*nh)
        yface = (nz-2*nh)*(nx-2*nh)
        xface = (ny-2*nh)*(nz-2*nh)

        faces = {'k': (nh,(nx-2*nh),(ny-2*nh)),
                 'j': (nh,(nz-2*nh),(nx-2*nh)),
                 'i': (nh,(ny-2*nh),(nz-2*nh))}

        self.sbuf = []
        self.rbuf = []
        self.reqs = []
        self.reqr = []
        self.nghbs = []
        self.direcs = []
        for direc in neighbours.keys():
            yourrank = neighbours[direc]
            if yourrank is None:
                pass
            else:
                shape = faces[direc[0]]
                sbuf = np.zeros(shape)
                rbuf = np.zeros(shape)
                reqs = comm.Send_init(sbuf, yourrank, myrank)
                reqr = comm.Recv_init(rbuf, yourrank, yourrank)

                self.sbuf += [sbuf]
                self.rbuf += [rbuf]
                self.reqs += [reqs]
                self.reqr += [reqr]
                self.nghbs += [yourrank]
                self.direcs += [direc]

    def fill(self, x):
        """ fill the halo of x"""

        reqs = self.reqs
        reqr = self.reqr

        MPI.Prequest.Startall(reqr)
        # halo to buffer
        # TODO: really copy the halo into the buffers
        # This operation is likely to be done in Fortran
        # We want to update the halo with only one pull of the array
        # from the memory
        for buf in self.sbuf:
            buf[:] = self.myrank

        MPI.Prequest.Startall(reqs)
        ierr = 9999
        try:
            ierr = MPI.Prequest.Waitall(reqr)
        except:
            print('IERR', ierr)
        # buffer to halo
        # TODO: fill the halo, really. Below is a dummy operation
        # this is the reverse operation (=>likely in Fortran)
        x[:] = 2.
        MPI.Prequest.Waitall(reqs)


if __name__ == '__main__':
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

    grid = {}
    grid['nh'] = nh
    grid['shape'] = (nz, ny, nx)
    grid['procs'] = procs

    msg = 'adjust the number of cores in mpirun -n XXX'
    assert nbproc == np.prod(procs), msg

    loc = rank2loc(myrank, procs)

    neighbours = get_neighbours(loc, procs)
    if myrank == 0:
        print('-'*40)
        print('neighbours:')
    comm.Barrier()
    print(myrank, loc, neighbours)
    comm.Barrier()

    halo = Halo_outerdim(comm, grid, neighbours)
#    halo = Halo(comm, grid, neighbours)
    print('halo is defined')
    comm.Barrier()

    x = np.zeros(grid['shape']) + myrank

    ntimes = 1000
    time = timeit.timeit(stmt='halo.fill(x)', number=ntimes,
                         globals={'halo': halo, 'x': x})
    comm.Barrier()
    for b, n, d in zip(halo.rbuf, halo.nghbs, halo.direcs):
        print('*', myrank, n, b[0, 0, 0], d)

    print(myrank, 'done')
    #print(myrank, loc, loc2rank(loc, procs))

    comm.Barrier()    
    if myrank == 0:        
        print('TIME = %.3e s per call' % (time/ntimes))
    
