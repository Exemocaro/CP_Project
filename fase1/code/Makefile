CC = gcc
SRC = src/
CFLAGS = -mavx -Wall -pg -O3 -march=native
# CFLAGS = -Wall -pg -O3 -march=native -mtune=native -flto -mavx -pthread
# -fopt-info-vec-missed -fopt-info-loop
# -fopt-info-optimized=output/report.txt
.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

clean:
	rm -f ./MD.exe gmon.out perf.data

run:
	./MD.exe < input/inputdata.txt
