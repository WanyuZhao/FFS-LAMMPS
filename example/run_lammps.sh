#!/bin/bash

#SBATCH -A 
#SBATCH -p 
#SBATCH --nodes=4
#SBATCH --ntasks=512
#SBATCH --time=3:00:00
#SBATCH --job-name test
#


mpirun -np 512 path/lmp_mpi -in lammps.input -screen none -ffs 128 ffs.input
