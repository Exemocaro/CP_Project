CC = gcc
MPICC = mpiCC
SRC = src/
CFLAGS = -Ofast -funroll-all-loops -ftree-vectorize -g

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe MDmpi.exe MDparmpi.exe

# Sequential version
MDseq.exe: $(SRC)/MDseq.cpp
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

# Parallel version with OpenMP
MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

# Parallel version with MPI
MDmpi.exe: $(SRC)/MDmpi.cpp
	$(MPICC) $(CFLAGS) $(SRC)MDmpi.cpp -lm -o MDmpi.exe

# Parallel version with MPI + OpenMP
MDparmpi.exe: $(SRC)/MDparmpi.cpp
	module load gcc/11.2.0;\
	$(MPICC) $(CFLAGS) $(SRC)MDparmpi.cpp -lm -fopenmp -o MDparmpi.exe

clean:
	rm -f ./MD*.exe

runseq: MDseq.exe
	sbatch scriptseq.sh

runpar: MDpar.exe
	sbatch scriptpar.sh

runmpi: MDmpi.exe
	sbatch scriptmpi.sh

runparmpi: MDparmpi.exe
	sbatch test_cluster2.sh 
