CC = gcc
SRC = src/
CFLAGS = -mavx -Wall -p -O3 -march=native
.DEFAULT_GOAL = all
# CFLAGS = -Wall -pg -O3 -march=native -mtune=native -flto -mavx -pthread
# -fopt-info-vec-missed -fopt-info-loop
# -fopt-info-optimized=output/report.txt

# Set the default number of threads to 8
#NUM_THREADS ?= 8

all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
#	module load gcc/11.2.0;
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
#	module load gcc/11.2.0;
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	rm -f ./MD*.exe gmon.out perf.data

runseq:
	./MDseq.exe < inputdata.txt

runpar:
	./MDpar.exe < inputdata.txt
