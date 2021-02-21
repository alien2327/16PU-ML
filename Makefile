CC = g++
CFLAGS = -IC:/tools/eigen -DEIGEN_STACK_ALLOCATION_LIMIT=0 -w
OBJS = main.o monitor.o genetic.o

all: main

main: monitor.o genetic.o main.o
	$(CC) $(OBJS) $(CFLAGS) -o main.exe

main.o: src/genetic.h src/monitor.h src/main.cpp
	$(CC) -c src/genetic.h src/monitor.h src/main.cpp $(CFLAGS)

monitor.o: src/monitor.h src/monitor.cpp
	$(CC) -c src/monitor.h src/monitor.cpp $(CFLAGS)

genetic.o: src/genetic.h src/genetic.cpp
	$(CC) -c src/genetic.h src/genetic.cpp $(CFLAGS)

.PHONY: clean

clean:
	del /f main.o monitor.o genetic.o