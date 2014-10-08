
objects = ballic.o tipsy.o
fortran_objects = icosahedron.o

CFLAGS ?= -O3

ballic: $(objects)
	cc -o ballic $(objects) $(fortran_objects) -lm

all: ballic

clean:
	rm $(objects)

