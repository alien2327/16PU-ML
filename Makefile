CC = g++
CFLAGS = -IC:/tools/eigen -DEIGEN_STACK_ALLOCATION_LIMIT=0 -w
OBJS = main.o monitor.o

all: main

main: monitor.o main.o
	$(CC) $(OBJS) $(CFLAGS) -o main.exe

main.o: src/monitor.h src/main.cpp
	$(CC) -c src/monitor.h src/main.cpp $(CFLAGS)

monitor.o: src/monitor.h src/monitor.cpp
	$(CC) -c src/monitor.h src/monitor.cpp $(CFLAGS)

.PHONY: clean

clean:
	del /f main.o monitor.o