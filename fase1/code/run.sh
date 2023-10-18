#!/bin/bash

sudo perf stat -e L1-dcache-load-misses,instructions,cycles -M cpi make run
cd output
diff cp_average.txt default_average.txt
diff cp_output.txt default_output.txt
cd ..