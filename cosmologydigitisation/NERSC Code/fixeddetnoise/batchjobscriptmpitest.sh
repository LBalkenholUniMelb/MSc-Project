#!/bin/bash
#SBATCH -N 2
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -J mpitest
#SBATCH --mail-user=lbalkenhol@student.unimelb.edu.au
#SBATCH --mail-type=ALL
#SBATCH -t 00:00:10

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
module load python/2.7-anaconda-4.4
srun -n 128 -c 1 --cpu_bind=threads python mpitest.py
