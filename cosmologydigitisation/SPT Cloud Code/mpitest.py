from mpi4py import MPI
comm = MPI.COMM_WORLD
processorrank = comm.Get_rank()

print(str(processorrank))