#!/bin/bash

sudo make clean
sudo make all
sudo perf stat -e cache-misses,cache-references,L1-dcache-load-misses,instructions,cycles -M cpi make runseq
cd output
diff cp_average.txt default_average.txt
diff cp_output.txt default_output.txt
cd ..