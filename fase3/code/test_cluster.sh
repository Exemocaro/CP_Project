#!/bin/bash
#SBATCH --time=1:00
#SBATCH --ntasks=3
#SBATCH --partition=cpar
mpirun -np 3 ./MDparmpi.exe < inputdata.txt
cd output

    diff ../cp_average.txt default_average.txt
    diff ../cp_output.txt default_output.txt   

    echo -n "average -> "
    diff -y --suppress-common-lines default_average.txt ../cp_average.txt | wc -l
    echo -n "output -> "
    diff -y --suppress-common-lines default_output.txt ../cp_output.txt | wc -l
cd ..
