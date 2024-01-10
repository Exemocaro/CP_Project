#!/bin/bash

make clean
#make all

make MDmpi.exe
make runmpi

cd output

    diff ../cp_average.txt default_average.txt
    diff ../cp_output.txt default_output.txt   

    echo -n "average -> "
    diff -y --suppress-common-lines default_average.txt ../cp_average.txt | wc -l
    echo -n "output -> "
    diff -y --suppress-common-lines default_output.txt ../cp_output.txt | wc -l
cd ..
