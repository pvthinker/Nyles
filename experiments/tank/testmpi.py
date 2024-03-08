from mpi4py import MPI


print(MPI.COMM_WORLD.Get_size())
MPI.Finalize()
