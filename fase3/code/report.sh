#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <function_name>"
    exit 1
fi

function_name="$1"

sudo perf record -g -- make runseq
sudo perf annotate --source "$function_name" > "output/report_${function_name}.txt"
