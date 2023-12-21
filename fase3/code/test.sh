#!/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --time=00:10:00
#SBATCH --partition=cpar
#SBATCH --exclusive

threads=(1 2 4 8 10 12 16 20 24 28 32 36 40)

for nthreads in "${threads[@]}"
do
	export OMP_NUM_THREADS=${nthreads}
	echo ${OMP_NUM_THREADS}
	time `./MDpar.exe <inputdata.txt >output/lixo`
done