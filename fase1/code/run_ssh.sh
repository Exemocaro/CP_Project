#!/bin/bash

srun --partition=cpar perf stat -e L1-dcache-load-misses,instructions,cycles make run
cd output
diff cp_average.txt default_average.txt
diff cp_output.txt default_output.txt
cd ..