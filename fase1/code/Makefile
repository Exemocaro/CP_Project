CC = gcc
SRC = src/
CFLAGS = -Wall -pg
# -Ofast -march=native -mtune=native -funroll-loops -flto -mavx -pthread

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

clean:
	rm ./MD.exe gmon.out

run:
	./MD.exe < input/inputdata.txt
