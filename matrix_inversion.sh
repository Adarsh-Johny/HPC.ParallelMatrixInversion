#!/bin/bash
#PBS -l select=2:ncpus=10:mem=2gb
#PBS -l walltime=0:01:00
#PBS -q short_cpuQ

# load the module
module load openmpi-4.0.4

# Actual running command
mpiexec ./HPC.ParallelMatrixInversion/main_program
