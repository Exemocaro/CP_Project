#!/bin/bash

sudo perf record -g -- make run
gprof MD.exe gmon.out > output/analysis.txt
gprof MD.exe gmon.out | gprof2dot | dot -Tpng -o images/out.png

