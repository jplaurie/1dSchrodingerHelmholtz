OBJS = main.o linear.o nonlinear.o read.o compute.o print.o forcing.o
#this is for linux machines
#CC = g++ -Wall -O3 -std=c++14 -fopenmp
#DEBUG = 
#CFLAGS = -c $(DEBUG)
#LFLAGS = $(DEBUG) -lfftw3 -lfftw3_threads -L/usr/local/lib/ -lopenblas -larmadillo -lm -lpthread

#The following is for mac os x
CC = g++-mp-13 -std=c++14 -Wall -O3 -fopenmp
DEBUG = #-ld_classic
CFLAGS = -c $(DEBUG) -I/opt/local/include
LFLAGS = $(DEBUG) -L/opt/local/lib -framework Accelerate -march=native -lfftw3 -lfftw3_threads -lopenblas -lm -lpthread -larmadillo



owt1d: $(OBJS)
	$(CC) -o owt1d $(OBJS) $(LFLAGS)

main.o: main.cpp const.h linear.h read.h print.h compute.h forcing.h
	$(CC) $(CFLAGS) main.cpp
	
nonlinear.o: nonlinear.cpp const.h
	$(CC) $(CFLAGS) nonlinear.cpp

linear.o: linear.cpp const.h
	$(CC) $(CFLAGS) linear.cpp

read.o: read.cpp const.h print.h
	$(CC) $(CFLAGS) read.cpp
    
compute.o: compute.cpp const.h nonlinear.h
	$(CC) $(CFLAGS) compute.cpp

print.o: print.cpp const.h
	$(CC) $(CFLAGS) print.cpp

forcing.o: forcing.cpp const.h
	$(CC) $(CFLAGS) forcing.cpp

clean:
	rm -rfv *.o ./output/ *~ owt1d

