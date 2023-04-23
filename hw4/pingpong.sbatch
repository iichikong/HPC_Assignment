#!/bin/bash
#SBATCH --job-name=pingpong
#SBATCH --output=pingpong.out
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=keigo.ando@nyu.edu

#Load the OpenMPI module.
module purge
module load intel/19.1.2 openmpi/intel/4.0.5

make

scontrol show hostname $SLURM_JOB_NODELIST >> pingpong.out
lscpu >> pingpong.out

# Run the pingpong program with mpi run
mpirun ./pingpong 0 1

make clean