#!/bin/bash

srun --partition=cpar perf stat -e L1-dcache-load-misses,instructions,cycles -M cpi make run