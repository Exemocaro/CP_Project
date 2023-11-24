

echo "###########################################################################"
echo "Running with 32 threads"
echo "###########################################################################"
export OMP_NUM_THREADS=32
sudo perf stat -e L1-dcache-load-misses -M cpi ./MDpar.exe < inputdata.txt
