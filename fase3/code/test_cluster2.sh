#!/bin/bash
#SBATCH --time=5:00
#SBATCH --ntasks=48
#SBATCH --partition=cpar

# Array of different numbers of processes to test
processes=(1 2 4 8 10 12 16 20 24 28 32 36 40 43 45 47)
#processes=(3 5)

# Loop through each number of processes
for nprocs in "${processes[@]}"
do
    echo "Running with ${nprocs} processes:"
    export OMP_NUM_THREADS=${processes}; time mpirun -np ${nprocs} ./MDparmpi.exe < inputdata.txt > output/output_${nprocs}.txt

    # cd output

    # diff ../cp_average.txt default_average.txt
    # diff ../cp_output.txt default_output.txt   

    # echo -n "average -> "
    # diff -y --suppress-common-lines default_average.txt ../cp_average.txt | wc -l
    # echo -n "output -> "
    # diff -y --suppress-common-lines default_output.txt ../cp_output.txt | wc -l
    # cd ..
    
done


