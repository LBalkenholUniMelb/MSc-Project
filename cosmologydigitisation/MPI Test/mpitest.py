from mpi4py import MPI
from numpy import *
from pickle import *

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

sendbuf = zeros((100, 100)) + rank
recvbuf = None

if rank == 0:
    recvbuf = empty((size, 100, 100))

comm.Gather(sendbuf, recvbuf, root = 0)

if rank == 0:
    conf = 0
    for prank, arr in recvbuf:
        if all(arr == prank):
            conf += 1
    confquot = float(conf)/float(size)

    dump(confquot, open("mpitestpickle.p", "wb"), -1)