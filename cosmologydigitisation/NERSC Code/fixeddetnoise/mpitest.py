from mpi4py import MPI

mpi_rank = MPI.COMM_WORLD.Get_rank()
mpi_size = MPI.COMM_WORLD.Get_size()

s = "rank: " + str(mpi_rank) + "::: size: " + str(mpi_size)
print(s)