.PHONY: all clean

all: main

main: main.o Mesh.o matrix.o
	g++ -o main main.o Mesh.o matrix.o

Mesh.o: Mesh.cpp Mesh.h
	g++ -c Mesh.cpp

matrix.o: matrix.cpp matrix.h
	g++ -c matrix.cpp

main.o: main.cpp
	g++ -c main.cpp

clean:
	rm *.o *.dat main

run:
	./main
