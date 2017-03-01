all: a.out

CC=g++
FDPS_PATH=$(HOME)/git/FDPS/src
CPPFLAGS=-I$(FDPS_PATH)
CPPFLAGS=-I/Users/watanabe/git/FDPS/src -std=c++11
CPPFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
LDFLAGS=-lmpi -lmpi_cxx

a.out: test.cpp
	$(CC) $(CPPFLAGS) $(LDFLAGS) $< -o $@

test: a.out
	mpirun -np 4 ./a.out
	gnuplot dump.plt

clean:
	rm -f a.out dump*.dat dump.png
