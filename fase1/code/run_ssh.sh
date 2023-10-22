#!/bin/bash

srun --partition=cpar perf stat -e L1-dcache-load-misses,instructions,cycles make run