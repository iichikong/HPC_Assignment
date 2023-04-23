#!/bin/bash
#SBATCH --job-name=jacobi
#SBATCH --output=jacobi.out
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=02:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=keigo.ando@nyu.edu

#Load the OpenMPI module.
module purge
module load intel/19.1.2 openmpi/intel/4.0.5

make

scontrol show hostname $SLURM_JOB_NODELIST >>jacobi.out
lscpu >>jacobi.out

# Run the pingpong program with mpi run
mpirun -np 16 ./jacobi 40 1000

make clean
