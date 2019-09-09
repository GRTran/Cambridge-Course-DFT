# Lines for GNU
CC=gcc
F90=gfortran
ELEC=-lhamiltonian
INCFLAGS=-I/usr/include/ 
LDFLAGS=-L. -Wl,-rpath=. -Bdynamic -llapack -lblas -lfftw3
TARGET=.

all: eigensolver

eigensolver: lib
	$(F90) -o eigensolver eigensolver.f90  $(ELEC)  -lm $(LDFLAGS)

lib: hamiltonian.f90
	$(F90) ${INCFLAGS}  -c -fPIC hamiltonian.f90 -lfftw3
	ar r libhamiltonian.a hamiltonian.o
	ranlib libhamiltonian.a
	$(F90)  -shared -Wl,-soname,libhamiltonian.so -o libhamiltonian.so  hamiltonian.o -lfftw3

clean:
	rm -f *.a *.o *.mod *.so
	rm eigensolver

