all: a.out

CC=g++
FDPS_PATH=$(HOME)/git/FDPS/src
CPPFLAGS=-I$(FDPS_PATH)
CPPFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
LDFLAGS=-lmpi -lmpi_cxx

a.out: test.cpp
	$(CC) $(CPPFLAGS) $(LDFLAGS) $< -o $@

test: a.out
	mpirun -np 2 ./a.out > test.dat
	gnuplot plot.plt

clean:
	rm -f a.out test.dat trac*.dat trajectory.png energy.png
