#!/bin/bash
#SBATCH --job-name=int_ring
#SBATCH --output=int_ring.out
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=keigo.ando@nyu.edu

#Load the OpenMPI module.
module purge
module load intel/19.1.2 openmpi/intel/4.0.5

make

scontrol show hostname $SLURM_JOB_NODELIST >>int_ring.out
lscpu >>int_ring.out

# Run the pingpong program with mpi run
mpirun ./int_ring 1000
mpirun ./int_ring_block 1000

make clean
