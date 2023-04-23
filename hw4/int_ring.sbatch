#!/bin/bash
#SBATCH --job-name=int_ring
#SBATCH --output=int_ring.out
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2GB
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=keigo.ando@nyu.edu

#Load the OpenMPI module.
module purge
module load intel/19.1.2 openmpi/intel/4.0.5

meke

scontrol show hostname $SLURM_JOB_NODELIST >>int_ring.out
lscpu >>int_ring.out

# Run the pingpong program with mpi run
mpirun ./int_ring 10
mpirun ./int_ring_block 10

make clean
