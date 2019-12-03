import numpy as np
import sys
from mpi4py import MPI

def get_myrank(procs):
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    msg = 'use mpirun -np %i python ' % np.prod(procs) + ' '.join(sys.argv)
    assert comm.Get_size() == np.prod(procs), msg

    return myrank

def abort():
    comm = MPI.COMM_WORLD
    comm.Abort()

def global_sum(localsum):
    return MPI.COMM_WORLD.allreduce(localsum, op=MPI.SUM)
