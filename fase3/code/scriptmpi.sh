#!/bin/bash
#SBATCH --time=5:00
#SBATCH --ntasks=40
#SBATCH --partition=cpar

time mpirun -np 40 ./MDmpi.exe < inputdata.txt
