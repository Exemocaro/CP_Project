CC = gcc
SRC = src/
CFLAGS = -Wall -pg -O3 -funroll-all-loops -ftree-vectorize -mavx

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
#	module load gcc/11.2.0;
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
#	module load gcc/11.2.0;
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	rm -f gmon.out perf.data ./MD*.exe

runseq: MDseq.exe
#	./MDseq.exe < inputdata.txt
	sudo perf stat -e cache-misses,cache-references,L1-dcache-load-misses,instructions,cycles ./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	export OMP_NUM_THREADS=4;\
	sudo perf stat -e cache-misses,cache-references,L1-dcache-load-misses,instructions,cycles ./MDpar.exe < inputdata.txt
#	./MDpar.exe < inputdata.txt
	
