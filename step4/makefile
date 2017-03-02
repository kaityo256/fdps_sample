all: a.out

CC=g++
FDPS_PATH=$(HOME)/git/FDPS/src
CPPFLAGS=-I$(FDPS_PATH)

a.out: test.cpp
	$(CC) $(CPPFLAGS) $(LDFLAGS) $< -o $@

test: a.out
	./a.out > test.dat
	gnuplot plot.plt

clean:
	rm -f a.out test.dat energy.png trajectory.png
