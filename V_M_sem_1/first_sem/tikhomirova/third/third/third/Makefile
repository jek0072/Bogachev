first: main.o otraz.o
	g++ otraz.o main.o -ffast-math -O3 -march=native -lpthread

otraz.o: otraz.cpp
	g++ -c otraz.cpp -ffast-math -O3 -march=native -lpthread
main.o: main.cpp
	g++ -c main.cpp -ffast-math -O3 -march=native -lpthread

