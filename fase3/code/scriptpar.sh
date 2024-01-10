#!/bin/bash
#SBATCH --time=5:00
#SBATCH --ntasks=40
#SBATCH --partition=cpar

export OMP_NUM_THREADS=40;\
	time ./MDpar.exe < inputdata.txt
