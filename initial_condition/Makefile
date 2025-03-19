OBJS = initial.o 
#this is for linux machines
CC = g++ -Wall -O3 -std=c++14 -fopenmp
#DEBUG = 
#CFLAGS = -c $(DEBUG) -I/usr/local/include/
#LFLAGS = $(DEBUG) -lfftw3 -lfftw3_threads -larmadillo -L/usr/local/lib/ -lm -lpthread

#this is for aston math cluster
#CC = g++ -Wall -O3 -std=c++14 -fopenmp
#CFLAGS = -c -I/apps/lib64/openblas/include -I/home/lauriej/libs/include -I/home/lauriej/libs/    usr/local/include
#LFLAGS = -lfftw3 -lfftw3_threads -L/appa/lib64/openbals/lib -L/home/lauriej/libs/lib -L/home/    lauriej/libs/usr/local/lib64 -lm -lpthread -lopenblas -llapack -larmadillo



#The following is for mac os x
CC = g++-mp-12 -Wall -O3 -std=c++14 -fopenmp
DEBUG = -ld_classic
CFLAGS = -c $(DEBUG) -I/opt/local/include
LFLAGS = $(DEBUG) -L/opt/local/lib -framework Accelerate -march=native -lfftw3 -lfftw3_threads -lopenblas -lm -lpthread -larmadillo



start: $(OBJS)
	$(CC) -o start $(OBJS) $(LFLAGS)

initial.o: initial.cpp
	$(CC) $(CFLAGS) initial.cpp
