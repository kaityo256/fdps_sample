all: a.out

CC=g++
FDPS_PATH=$(HOME)/git/FDPS/src
CPPFLAGS=-I$(FDPS_PATH)
CPPFLAGS += -std=c++11
CPPFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
LDFLAGS=-lmpi -lmpi_cxx

a.out: test.cpp
	$(CC) $(CPPFLAGS) $(LDFLAGS) $< -o $@

test: a.out
	mpirun -np 4 ./a.out
	gnuplot dump.plt
	convert plot*.png plot.gif

clean:
	rm -f a.out dump*.dat *.png plot.gif
