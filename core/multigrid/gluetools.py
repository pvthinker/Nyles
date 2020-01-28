"""
Tools to glue/split arrays and matrices together

fine  ===Restrict==> dummy ===Glue===> coarse
coarse ====Split===> dummy ===Interpol===> fine

where Restrict and Interpol are matrix-vector multiplications,
handled by 'intergrids'

dummy grid is defined only if the gluing operation is required
dummy grid has the same neighbours than fine

"""
import numpy as np
from mpi4py import MPI
import grids as gr
import topology as topo
from scipy import sparse
import pickle
import os
import itertools
import mpitools


class Gluegrids(object):
    def __init__(self, fine, coarse):
        self.set_localcomm(fine, coarse)
        self.set_dummy(fine, coarse)
        self.fine = fine
        self.coarse = coarse

    def set_localcomm(self, fine, coarse):
        """ define a local communicator for the gluing

        among this communicator, all cores will do exactly them
        same work from this level to the coarsest.

        we mostly use the Allgatherv instruction, except for the
        matrix gluing where it turns out to be easier to implement
        a core-to-core communication.
        """
        comm = MPI.COMM_WORLD
        myrank = comm.Get_rank()
        procs0 = [i*p for i, p in zip(coarse.incr, coarse.procs)]
        mat = topo.get_mypartners(procs0, fine.incr, coarse.incr, myrank)
        self.matshape = np.shape(mat)
        partners = mat.flat[:]
        master = partners[0]

        self.myrank = myrank
        self.mat = mat

        # http://mpi.deino.net/mpi_functions/MPI_Comm_split.html
        self.localcomm = MPI.COMM_WORLD.Split(master, myrank)

    def set_dummy(self, fine, coarse):
        matshape = self.matshape
        self.dummy = gr.Dummygrid(fine, coarse)
        # exchange domainindices
        msg = np.array(self.dummy.domainindices, dtype=int)
        length = len(msg)
        rbuff = np.zeros(matshape+(length,), dtype=int)
        siz = np.ones(np.prod(matshape))*length
        off = np.zeros(np.prod(matshape))
        off[1:] = np.cumsum(siz)[:-1]
        self.localcomm.Allgatherv(msg, [rbuff, siz, off, MPI.INT64_T])
        self.domainindices = rbuff.reshape((np.prod(matshape), length))

        # exchange sizes
        msg = np.array(self.dummy.size, dtype=int)
        length = len(msg)
        rbuff = np.zeros(matshape+(length,), dtype=int)
        siz = np.ones(np.prod(matshape))*length
        off = np.zeros(np.prod(matshape))
        off[1:] = np.cumsum(siz)[:-1]
        self.localcomm.Allgatherv(msg, [rbuff, siz, off, MPI.INT64_T])
        self.sizes = rbuff.reshape((np.prod(matshape), length))

        # pre-allocate buffers for glue_array
        dummy = self.dummy
        matshape = self.matshape
        # exchange vectors (glue dummy vectors together)
        k0, k1, j0, j1, i0, i1 = dummy.domainindices
        msg = dummy.toarray("x")[k0:k1, j0:j1, i0:i1].ravel()
        length = len(msg)
        self.rbuff = np.zeros(matshape+(length,))

    def set_dummy_Amatrix(self):
        """ for testing/debugging purposes only
        in reality dummy.A is defined in mg.py """
        dummy = self.dummy
        N = dummy.N
        k0, k1, j0, j1, i0, i1 = dummy.domainindices
        dummy.x[:] = 0.
        dummy.toarray("x")[k0:k1, j0:j1, i0:i1] = 1.
        dummy.A = sparse.spdiags(dummy.x, 0, N, N)

    def glue_array(self, which):
        """ glue arrays dummy.x into coarse.which

        dummy.x needs to be assigned before this stage
        it is assigned by intergrids.py by applying the
        restriction matrix on fine.x

        """
        assert which in "xrb"

        coarse = self.coarse
        dummy = self.dummy
        matshape = self.matshape
        # exchange vectors (glue dummy vectors together)
        k0, k1, j0, j1, i0, i1 = dummy.domainindices
        msg = dummy.toarray("x")[k0:k1, j0:j1, i0:i1].ravel()
        length = len(msg)
        rbuff = self.rbuff  # pre-allocated
        siz = np.ones(np.prod(matshape))*length
        off = np.zeros(np.prod(matshape))
        off[1:] = np.cumsum(siz)[:-1]
        self.localcomm.Allgatherv(msg, [rbuff, siz, off, MPI.DOUBLE])
        x = rbuff.reshape((np.prod(matshape),) + tuple(dummy.shape))

        l = 0
        nz, ny, nx = dummy.shape
        coarseb = coarse.toarray(which)
        k0, k1, j0, j1, i0, i1 = coarse.domainindices
        for k in range(matshape[0]):
            ka = k0+k*nz
            kb = ka+nz
            for j in range(matshape[1]):
                ja = j0+j*ny
                jb = ja+ny
                for i in range(matshape[2]):
                    ia = i0+i*nx
                    ib = ia+nx
                    coarseb[ka:kb, ja:jb, ia:ib] = x[l][:]
                    l += 1
    # the halofill is done in intergrids
    # coarse.halofill(which)

    def split_array(self, which):
        """ split array coarse.which into dummy.x

        """
        assert which in "xrb"
        x = self.dummy.toarray("x")
        xcoarse = self.coarse.toarray(which)
        #
        k0, k1, j0, j1, i0, i1 = self.coarse.domainindices
        k, j, i = np.where(self.myrank == self.mat)
        k, j, i = int(k), int(j), int(i)
        nz, ny, nx = self.dummy.shape
        ka, kb = k0+k*nz, k0+k*nz+nz
        ja, jb = j0+j*ny, j0+j*ny+ny
        ia, ib = i0+i*nx, i0+i*nx+nx
        #
        k0, k1, j0, j1, i0, i1 = self.dummy.domainindices
        x[k0:k1, j0:j1, i0:i1] = xcoarse[ka:kb, ja:jb, ia:ib]

    def glue_matrix(self):
        """ glue dummy.A into coarse.A

        dummy.A is assigned outside of this function
        """
        coarse = self.coarse
        dummy = self.dummy
        matshape = self.matshape
        myrank = self.myrank
        mat = self.mat

        comm = MPI.COMM_WORLD
        # exchange matrices

        N = coarse.N
        row = np.zeros((N*27,))
        col = np.zeros((N*27,))
        data = np.zeros((N*27,))
        counter = 0
        k0, k1, j0, j1, i0, i1 = coarse.domainindices
        nz, ny, nx = dummy.shape
        ok = True
        rks = range(matshape[0])
        rjs = range(matshape[1])
        ris = range(matshape[2])
        master = mat[0, 0, 0]
        for k, j, i in itertools.product(rks, rjs, ris):
            ka = k0+k*nz
            ja = j0+j*ny
            ia = i0+i*nx
            # core to core communication
            sender = mat[k, j, i]
            if myrank == sender:
                obj = {"A": dummy.A.tocoo(),
                       "size": dummy.size,
                       "domainindices": dummy.domainindices}
                for dest in mat.ravel():
                    if dest != myrank:
                        # in python MPI can send any dictionary
                        # this makes the matrix exchange almost easy
                        # the same implementation in Fortran
                        # would be quite a nightmare
                        comm.send(obj, dest=dest, tag=myrank)
            else:
                obj = comm.recv(source=sender, tag=sender)
            comm.Barrier()

            A = obj["A"]
            size = obj["size"]
            domainindices = obj["domainindices"]

            n = len(A.row)

            # row
            kk, jj, ii = index_to_triplet(A.row, size, domainindices)

            kk += ka-k0
            jj += ja-j0
            ii += ia-i0
            idx = triplet_to_index(
                (kk, jj, ii), coarse.size, coarse.domainindices)
            row[counter:counter+n] = idx
            if (max(idx) >= N):
                print("myrank=%i, sender=%i / row" % (myrank, sender))

                ok = False

            # col
            kk, jj, ii = index_to_triplet(A.col, size, domainindices)
            kk += ka-k0
            jj += ja-j0
            ii += ia-i0
            idx = triplet_to_index(
                (kk, jj, ii), coarse.size, coarse.domainindices)

            col[counter:counter+n] = idx
            if (max(idx) >= N):
                print("myrank=%i, sender=%i / col" % (myrank, sender))
                ok = False

            # data
            data[counter:counter+n] = A.data

            counter += n

        MPI.COMM_WORLD.Barrier()
        if ok:
            coarseA = sparse.coo_matrix((data, (row, col)), shape=(N, N))

        else:
            print("matrix not defined for rank %i" % myrank)
            print("subdomains were", mat)
            raise ValueError

        return coarseA


def triplet_to_index(triplet, size, domainindices):
    """ 
    (0,0,0) is the bottom left point in the *interior*
    (0,0,-1) is in the halo
    (0,0,nx) is in the halo
    """
    k, j, i = triplet
    k0, k1, j0, j1, i0, i1 = domainindices
    nnz, nny, nnx = size
    index = (i+i0) + (j+j0)*nnx + (k+k0)*nnx*nny
    return index


def index_to_triplet(index, size, domainindices):
    k0, k1, j0, j1, i0, i1 = domainindices
    nnz, nny, nnx = size
    k = index // (nnx*nny)
    j = (index-k*(nnx*nny)) // nnx
    i = (index-k*(nnx*nny)-j*nnx)
    k -= k0
    j -= j0
    i -= i0
    return (k, j, i)
