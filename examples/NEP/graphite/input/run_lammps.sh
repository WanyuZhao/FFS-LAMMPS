#!/bin/bash
#PBS -A 
#PBS -q 
#PBS -l select=4:ncpus=128:mpiprocs=128
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N 3900-15



aprun -n 512 path/src/lmp_mpi -in lammps.input -screen none -ffs 128 ffs.input 

