OBJS = main.o linear.o nonlinear.o read.o energy.o intensity.o spectrum.o initialise.o flux.o
#this is for linux machines
CC = g++ -Wall -O3 -std=c++14 -fopenmp
DEBUG = 
CFLAGS = -c $(DEBUG)
LFLAGS = $(DEBUG) -lfftw3 -lfftw3_threads -L/usr/local/lib/ -lopenblas -larmadillo -lm -lpthread

#The following is for mac os x
#CC = g++-mp-11 -Wall -O3 -fopenmp
#DEBUG = 
#CFLAGS = -c $(DEBUG) -I/opt/local/include
#LFLAGS = $(DEBUG) -L/opt/local/lib -framework Accelerate -march=native -lfftw3 -lfftw3_threads -lopenblas -lm -lpthread -larmadillo



owt1d-long: $(OBJS)
	$(CC) -o owt1d-long $(OBJS) $(LFLAGS)

main.o: main.cpp const.h functions.h
	$(CC) $(CFLAGS) main.cpp
	
nonlinear.o: nonlinear.cpp const.h functions.h
	$(CC) $(CFLAGS) nonlinear.cpp

linear.o: linear.cpp const.h functions.h
	$(CC) $(CFLAGS) linear.cpp

read.o: read.cpp const.h functions.h
	$(CC) $(CFLAGS) read.cpp
    
energy.o: energy.cpp const.h functions.h
	$(CC) $(CFLAGS) energy.cpp

intensity.o: intensity.cpp const.h functions.h
	$(CC) $(CFLAGS) intensity.cpp

spectrum.o: spectrum.cpp const.h functions.h
	$(CC) $(CFLAGS) spectrum.cpp

#Forcing.o: Forcing.cpp Const.h Functions.h
#	$(CC) $(CFLAGS) Forcing.cpp

initialise.o: initialise.cpp const.h functions.h
	$(CC) $(CFLAGS) initialise.cpp

#Filter.o: Filter.cpp Const.h Functions.h
#	$(CC) $(CFLAGS) Filter.cpp

flux.o: flux.cpp const.h
	$(CC) $(CFLAGS) flux.cpp

clean:
	rm -rfv *.o ./output/ *~ owt1d-long

