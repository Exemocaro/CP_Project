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
#	module load gcc/11.2.0;\.
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

# Parallel version with MPI
MDmpi.exe: $(SRC)/MDmpi.cpp
	$(MPICC) $(CFLAGS) $(SRC)MDmpi.cpp -lm -o MDmpi.exe

# Parallel version with MPI + OpenMP
MDparmpi.exe: $(SRC)/MDparmpi.cpp
#   module load gcc/11.2.0;\.
	$(MPICC) $(CFLAGS) $(SRC)MDparmpi.cpp -lm -fopenmp -o MDparmpi.exe

clean:
	rm -f ./MD*.exe

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	export OMP_NUM_THREADS=40;\
	./MDpar.exe < inputdata.txt

runmpi: MDmpi.exe
	mpirun -np 4 ./MDmpi.exe < inputdata.txt

runparmpi: MDparmpi.exe
	export OMP_NUM_THREADS=40;\
	mpirun -np 4 ./MDparmpi.exe < inputdata.txt
