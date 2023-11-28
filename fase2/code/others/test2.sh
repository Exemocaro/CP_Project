#!/bin/bash
#SBATCH --time=2:00
#SBATCH --ntasks=40
#SBATCH --partition=cpar

echo "###########################################################################"
echo "Running with 32 threads"
echo "###########################################################################"
export OMP_NUM_THREADS=32
perf stat -e L1-dcache-load-misses -M cpi ./MDpar.exe < inputdata.txt
