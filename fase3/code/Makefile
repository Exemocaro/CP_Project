CC = gcc
MPICC = mpiCC
SRC = src/
CFLAGS = -Ofast -funroll-all-loops -ftree-vectorize -g

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe MDparmpi.exe

# Sequential version
MDseq.exe: $(SRC)/MDseq.cpp
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

# Parallel version with OpenMP
MDpar.exe: $(SRC)/MDpar.cpp
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

# Parallel version with MPI
MDparmpi.exe: $(SRC)/MDparmpi.cpp
	$(MPICC) $(CFLAGS) $(SRC)MDparmpi.cpp -lm -fopenmp -o MDparmpi.exe

clean:
	rm -f ./MD*.exe

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	export OMP_NUM_THREADS=40;\
	./MDpar.exe < inputdata.txt

runparmpi: MDparmpi.exe
	export OMP_NUM_THREADS=40;\
	mpirun -np 3 ./MDparmpi.exe < inputdata.txt
