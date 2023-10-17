CC = gcc
SRC = src/
CFLAGS = -ftree-vectorize -msse2 -Wall -g -O3 
# CFLAGS = -Wall -pg -O3 -march=native -mtune=native -flto -mavx -pthread
# -fopt-info-vec-missed
.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

clean:
	rm -f ./MD.exe gmon.out

run:
	./MD.exe < input/inputdata.txt
