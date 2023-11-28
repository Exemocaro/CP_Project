#!/bin/bash

sudo make clean
sudo make all
#sudo make runseq
sudo make runpar

cd output
    echo -n "average -> "
    diff -y --suppress-common-lines default_average.txt cp_average.txt | wc -l
    echo -n "output -> "
    diff -y --suppress-common-lines default_output.txt cp_output.txt | wc -l
cd ..