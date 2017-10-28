CC=g++
CFLAGS=-c -Wall -I. -std=c++11

all: model

model: main.o model.o potential.o
	$(CC) main.o model.o potential.o -o model

model.o: model.cpp
	$(CC) $(CFLAGS) model.cpp

potential.o: potential.cpp
	$(CC) $(CFLAGS) potential.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

clean:
	rm -rf main.o model.o potential.o model
