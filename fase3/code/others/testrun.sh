#!/bin/bash

# SLURM directives
#SBATCH --ntasks=40
#SBATCH --time=00:10:00
#SBATCH --partition=cpar
#SBATCH --exclusive

# Clean and build the project
sudo make clean
sudo make all

# Run the application with different thread numbers
threads=(1 2 4 20 40)
for nthreads in "${threads[@]}"
do
    export OMP_NUM_THREADS=${nthreads}
    echo "Running with OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    time ./MDpar.exe <inputdata.txt >lixo
done

# Performance profiling
sudo perf stat -e cache-misses,cache-references,L1-dcache-load-misses,instructions,cycles -M cpi make runpar

# Comparing output files
cd output
diff cp_average.txt default_average.txt
diff cp_output.txt default_output.txt
cd ..
